from collections import defaultdict
from operator import itemgetter
import warnings

import pandas as pd
import pybedtools as pbt
import numpy as np
from scipy.stats import rankdata

GTF_EXTENSIONS = (
    ".gtf",
    ".gtf.gz",
)


def _compute_union_exon_length(starts_ends):
    """Compute gene length based on the union of its exons.
    Before summing up exon lengths, one needs to merge all overlapping
    exons. To do this, ``pybedtools.merge()`` could be used, but simple
    hand-built merge works much faster.
    This function works only if start_ends is a list of lists. It will
    not work with list of tuples.
    """
    starts_ends.sort(key=itemgetter(0))
    merged = [starts_ends[0]]
    for current in starts_ends:
        previous = merged[-1]
        if current[0] <= previous[1]:
            previous[1] = max(previous[1], current[1])
        else:
            merged.append(current)

    return sum(end - start for (start, end) in merged)


def union_exon_lengths(annotation, gene_id_attr="gene_id"):
    """Compute gene lengths based on union exon model of genome annotation.
    Group exon start & end coordinates by gene ID (level 1) and chrom &
    strand (level 2). Then perfrom merge and length calculation for each
    chrom & strand separately. The latter is needed since
    ``gene_id_attr`` is not unique in some annotations (e.g. RefSeq).
    This function is confirmed to produce identical output as
    featureCounts output (column "Length") for the following species /
    annotations:
      - Homo sapiens:
        - UCSC hg19
        - UCSC hg38
        - ENSEMBL 92
        - ENSEMBL 100
      - Mus musculus:
        - UCSC mm10
        - ENSEMBL 92
        - ENSEMBL 100
      - Rattus norvegicus:
        - UCSC rn6
        - ENSEMBL 92
        - ENSEMBL 100
      - Macaca mulatta:
        - ENSEMBL 97
        - ENSEMBL 100
    """
    if not annotation.endswith(GTF_EXTENSIONS):
        raise ValueError(f"Input file ({annotation}) should be in GTF format")

    data = defaultdict(lambda: defaultdict(list))
    for segment in pbt.BedTool(annotation):
        if segment[2] != "exon":
            continue
        if gene_id_attr not in segment.attrs:
            raise ValueError(
                "Gene ID attribute is missing in the segment {segment[:]}. Please "
                "supply correct gene ID attribute with --gene-id-attr parameter."
            )

        data[segment.attrs[gene_id_attr]][(segment.chrom, segment.strand)].append([segment.start, segment.end])

    gene_lengths = defaultdict(int)
    for gene_id, by_chrom_strand in data.items():
        for starts_ends in by_chrom_strand.values():
            gene_lengths[gene_id] += _compute_union_exon_length(starts_ends)

    df = pd.DataFrame.from_dict(gene_lengths, orient="index", columns=["GENE_LENGTHS"])
    df.index.names = ["FEATURE_ID"]
    return df

def _tpm_ndarray(X, y):
    """Normalize Numpy ndarray expression counts to TPM.
    :type X: Numpy ndarray
    :type y: Numpy ndarray
    """
    assert isinstance(X, np.ndarray)
    assert isinstance(y, np.ndarray)
    assert X.shape[0] == y.shape[0]
    assert y.shape[1] == 1
    assert np.min(X) >= 0.0  # Gene counts must be non-negative

    A = X / y
    sumA = A.sum(axis=0)

    with np.errstate(invalid="ignore"):  # Ignore warnings of division by 0
        TPM = 1e6 * A / sumA

        # Samples with zeros for all genes get nan but should be 0.0
        np.nan_to_num(TPM, copy=False)

    return TPM


def tpm(X, y):
    """Normalize expression counts to Transcript per kilobase million (TPM).
    A = readsMappedToGene / geneLength
    TPM = A / SUM(A) * 1e6
    :type X: 2-D array_like
    :type y: 1-D array_like
    """
    if isinstance(X, pd.DataFrame) and isinstance(y, pd.DataFrame):
        common_genes = X.index.intersection(y.index)
        ncommon = len(common_genes)
        if ncommon != X.shape[0] or ncommon != y.shape[0]:
            warnings.warn(
                f"Geneset mismatch between expressions ({X.shape[0]}) and gene "
                f"lengths ({y.shape[0]}). Using intersection ({ncommon})...",
                RuntimeWarning,
                stacklevel=2,
            )

        X = X.loc[common_genes]
        y = y.loc[common_genes]

        X_ = np.asarray(X, dtype=np.float64)
        y_ = np.asarray(y, dtype=np.float64)
        TPM = _tpm_ndarray(X_, y_)

        return pd.DataFrame(TPM, index=X.index, columns=X.columns)
    else:
        X = np.asarray(X, dtype=np.float64)
        y = np.asarray(y, dtype=np.float64)
        return _tpm_ndarray(X, y)


def _fpkm_ndarray(X, y):
    """Normalize Numpy ndarray expression counts to FPKM.
    :type X: Numpy ndarray
    :type y: Numpy ndarray
    """
    assert isinstance(X, np.ndarray)
    assert isinstance(y, np.ndarray)
    assert X.shape[0] == y.shape[0]
    assert y.shape[1] == 1
    assert np.min(X) >= 0.0  # Gene counts must be non-negative

    total_sample_reads = X.sum(axis=0) / 10**6

    with np.errstate(invalid="ignore"):  # Ignore warnings of division by 0
        rpm = X / total_sample_reads
        # Samples with zeros for all genes get nan but should be 0.0
        np.nan_to_num(rpm, copy=False)

    fpkm = rpm / y * 1000

    return fpkm


def fpkm(X, y):
    """Normalize expression counts to Fragments per kilobase million (FPKM).
    RPM = expressionCount / sumReadsInSample * 1e6
    FPKM = RPM / geneLengthKb
    :type X: 2-D array_like
    :type y: 1-D array_like
    """
    if isinstance(X, pd.DataFrame) and isinstance(y, pd.DataFrame):
        common_genes = X.index.intersection(y.index)
        ncommon = len(common_genes)
        if ncommon != X.shape[0] or ncommon != y.shape[0]:
            warnings.warn(
                f"Geneset mismatch between expressions ({X.shape[0]}) and gene "
                f"lengths ({y.shape[0]}). Using intersection ({ncommon})...",
                RuntimeWarning,
                stacklevel=2,
            )

        X = X.loc[common_genes]
        y = y.loc[common_genes]

        X_ = np.asarray(X, dtype=np.float64)
        y_ = np.asarray(y, dtype=np.float64)
        FPKM = _fpkm_ndarray(X_, y_)

        return pd.DataFrame(FPKM, index=X.index, columns=X.columns)
    else:
        X = np.asarray(X, dtype=np.float64)
        y = np.asarray(y, dtype=np.float64)
        return _fpkm_ndarray(X, y)


def cpm(X):
    """Normalize expression counts to Counts per million (CPM).
    CPM = readsMappedToGene / totalNumReads * 1e6
    :type X: 2-D array_like
    """
    if isinstance(X, pd.DataFrame):
        assert X.to_numpy().min() >= 0.0  # Gene counts must be non-negative
    else:
        # Cast non-Pandas array_like objects to Numpy
        X = np.asarray(X, dtype=np.float64)
        assert np.min(X) >= 0.0  # Gene counts must be non-negative

    sumX = X.sum(axis=0)

    with np.errstate(invalid="ignore"):  # Ignore warnings of division by 0
        CPM = 1e6 * X / sumX

        # Samples with zeros for all genes get nan but should be 0.0
        np.nan_to_num(CPM, copy=False)

    return CPM


def quantile(X):
    """Quantile normalize gene expression to average distribution.
    The procedure is implemented as described on Wikipedia_ and runs on columns and rows:
    * Rearrange the column values so each column is in order from lowest to highest value
    * Find the mean for each row to determine the average distribution of expression values
    * For each column in original data determine a rank from lowest to highest
    * Take the ranking order and substitute in new values from the average distribution
    .. _Wikipedia: https://en.wikipedia.org/wiki/Quantile_normalization
    """
    # Cast array_like objects to Numpy
    X_ = np.asarray(X, dtype=np.float64)
    assert np.min(X_) >= 0.0  # Gene expression must be non-negative

    average_expression_distribution = np.mean(np.sort(X_, axis=0), axis=1)

    rank_avg = rankdata(X_, method="average", axis=0) - 1
    rank_floor = rank_avg.astype(int)
    rank_ceil = np.ceil(rank_avg).astype(int)

    X_floor = average_expression_distribution.take(rank_floor)
    X_ceil = average_expression_distribution.take(rank_ceil)
    X_ = (X_floor + X_ceil) / 2

    if isinstance(X, pd.DataFrame):
        return pd.DataFrame(X_, index=X.index, columns=X.columns)
    else:
        return X_