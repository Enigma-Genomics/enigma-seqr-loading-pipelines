from collections.abc import Callable

import hail as hl

from loading_pipeline.lib.annotations.enums import (
    BIOTYPES,
    FIVEUTR_CONSEQUENCES,
    LOF_FILTERS,
    TRANSCRIPT_CONSEQUENCE_TERMS,
)

BIOTYPE_LOOKUP = hl.dict(hl.enumerate(BIOTYPES, index_first=False))
EXTENDED_INTRONIC_SPLICE_REGION_VARIANT = 'extended_intronic_splice_region_variant'
FIVEUTR_CONSEQUENCES_LOOKUP = hl.dict(
    hl.enumerate(FIVEUTR_CONSEQUENCES, index_first=False).extend(
        [(hl.missing(hl.tstr), hl.missing(hl.tint32))],
    ),
)
LOF_FILTERS_LOOKUP = hl.dict(hl.enumerate(LOF_FILTERS, index_first=False))
MANE_SELECT_ANNOTATIONS = [
    'mane_select',
    'mane_plus_clinical',
]
NAGNAG_SITE = 'NAGNAG_SITE'
OMIT_TRANSCRIPT_CONSEQUENCE_TERMS = hl.set(
    [
        'upstream_gene_variant',
        'downstream_gene_variant',
    ],
)
PROTEIN_CODING_ID = BIOTYPE_LOOKUP['protein_coding']
SELECTED_ANNOTATIONS = [
    'amino_acids',
    'canonical',
    'codons',
    'gene_id',
    'hgvsc',
    'hgvsp',
    'transcript_id',
]
TRANSCRIPT_CONSEQUENCE_TERMS_LOOKUP = hl.dict(
    hl.enumerate(TRANSCRIPT_CONSEQUENCE_TERMS, index_first=False),
)


def _lof_filter_ids(c: hl.StructExpression) -> hl.ArrayNumericExpression:
    return hl.or_missing(
        (c.lof == 'LC') & hl.is_defined(c.lof_filter),
        c.lof_filter.split('&|,').map(lambda f: LOF_FILTERS_LOOKUP[f]),
    )


def _consequence_term_ids(c: hl.StructExpression) -> hl.ArrayNumericExpression:
    return c.consequence_terms.filter(
        lambda t: ~OMIT_TRANSCRIPT_CONSEQUENCE_TERMS.contains(t),
    ).map(lambda t: TRANSCRIPT_CONSEQUENCE_TERMS_LOOKUP[t])


def _legacy_utr_effects(c):
    # Retain the information old records actually contain, including every entry.
    return hl.or_else(hl.or_missing(hl.is_defined(c.fiveutr_consequence), hl.array([
        hl.struct(
            effect=hl.missing(hl.tstr),
            consequence=hl.or_missing(hl.is_defined(FIVEUTR_CONSEQUENCES_LOOKUP.get(c.fiveutr_consequence)), c.fiveutr_consequence),
            annotations=hl.sorted(c.fiveutr_annotation.items(), key=lambda entry: entry[0]).map(
                lambda entry: hl.struct(index=entry[0], annotation=entry[1]),
            ),
        ),
    ])), hl.empty_array(hl.tstruct(
        effect=hl.tstr, consequence=hl.tstr,
        annotations=hl.tarray(hl.tstruct(index=hl.tstr, annotation=c.fiveutr_annotation.dtype.value_type)),
    )))


def _all_utr_effects(c):
    # Sort classes and annotation entries for serialization, never severity.
    if 'fiveutr_all_effects' not in c.dtype.fields:
        return _legacy_utr_effects(c)
    effects = hl.sorted(
        hl.or_else(c.fiveutr_all_effects.items(), hl.empty_array(
            hl.ttuple(hl.tstr, c.fiveutr_all_effects.dtype.value_type),
        )), key=lambda item: item[0],
    ).map(lambda item: hl.struct(
        effect=item[0],
        consequence=hl.or_missing(hl.is_defined(FIVEUTR_CONSEQUENCES_LOOKUP.get(item[1].consequence)), item[1].consequence),
        annotations=hl.sorted(item[1].annotation.items(), key=lambda entry: entry[0]).map(
            lambda entry: hl.struct(index=entry[0], annotation=entry[1]),
        ),
    ))
    return hl.if_else(hl.is_defined(c.fiveutr_all_effects), effects, _legacy_utr_effects(c))


def vep_110_transcript_consequences_select(
    gencode_ensembl_to_refseq_id_mapping: hl.tdict(hl.tstr, hl.tstr),
) -> hl.StructExpression:
    return lambda c: c.select(
        *SELECTED_ANNOTATIONS,
        *MANE_SELECT_ANNOTATIONS,
        biotype_id=BIOTYPE_LOOKUP[c.biotype],
        consequence_term_ids=_consequence_term_ids(c),
        exon=hl.bind(
            lambda split: hl.or_missing(
                hl.is_defined(split),
                hl.Struct(index=split[0], total=split[1]),
            ),
            c.exon.split('/').map(hl.parse_int32),
        ),
        intron=hl.bind(
            lambda split: hl.or_missing(
                hl.is_defined(split),
                hl.Struct(index=split[0], total=split[1]),
            ),
            c.intron.split('/').map(hl.parse_int32),
        ),
        refseq_transcript_id=gencode_ensembl_to_refseq_id_mapping.get(c.transcript_id),
        alphamissense=hl.struct(
            pathogenicity=c.am_pathogenicity,
        ),
        loftee=hl.struct(
            is_lof_nagnag=c.lof_flags == NAGNAG_SITE,
            lof_filter_ids=_lof_filter_ids(c),
        ),
        spliceregion=hl.struct(
            extended_intronic_splice_region_variant=(
                hl.is_defined(c.spliceregion)
                & c.spliceregion.contains(
                    EXTENDED_INTRONIC_SPLICE_REGION_VARIANT,
                )
            ),
        ),
        utrannotator=hl.struct(
            fiveutr_consequences=_all_utr_effects(c).map(lambda effect: effect.consequence),
            fiveutr_effects_json=hl.json(_all_utr_effects(c)),
            existing_inframe_oorfs=c.existing_inframe_oorfs,
            existing_outofframe_oorfs=c.existing_outofframe_oorfs,
            existing_uorfs=c.existing_uorfs,
            fiveutr_consequence_id=FIVEUTR_CONSEQUENCES_LOOKUP[c.fiveutr_consequence],
            # Annotation documentation here:
            # https://github.com/ImperialCardioGenetics/UTRannotator?tab=readme-ov-file#the-detailed-annotation-for-each-consequence
            # NB:
            fiveutr_annotation=c.fiveutr_annotation['1'].annotate(
                AltStopDistanceToCDS=hl.parse_int32(
                    c.fiveutr_annotation['1'].AltStopDistanceToCDS,
                ),
                CapDistanceToStart=hl.parse_int32(
                    c.fiveutr_annotation['1'].CapDistanceToStart,
                ),
                DistanceToCDS=hl.parse_int32(
                    c.fiveutr_annotation['1'].DistanceToCDS,
                ),
                DistanceToStop=hl.parse_int32(
                    c.fiveutr_annotation['1'].DistanceToStop,
                ),
                Evidence=hl.or_missing(
                    # Just in case a weird value ("NA" or anything else) propagates
                    (
                        (c.fiveutr_annotation['1'].Evidence == 'True')
                        | (c.fiveutr_annotation['1'].Evidence == 'False')
                    ),
                    hl.bool(c.fiveutr_annotation['1'].Evidence),
                ),
                StartDistanceToCDS=hl.parse_int32(
                    c.fiveutr_annotation['1'].StartDistanceToCDS,
                ),
                newSTOPDistanceToCDS=hl.parse_int32(
                    c.fiveutr_annotation['1'].newSTOPDistanceToCDS,
                ),
                alt_type_length=hl.parse_int32(
                    c.fiveutr_annotation['1'].alt_type_length,
                ),
                ref_StartDistanceToCDS=hl.parse_int32(
                    c.fiveutr_annotation['1'].ref_StartDistanceToCDS,
                ),
                ref_type_length=hl.parse_int32(
                    c.fiveutr_annotation['1'].ref_type_length,
                ),
            ),
        ),
    )


def vep_85_transcript_consequences_select(
    c: hl.StructExpression,
) -> hl.StructExpression:
    return c.select(
        *SELECTED_ANNOTATIONS,
        biotype_id=BIOTYPE_LOOKUP[c.biotype],
        consequence_term_ids=_consequence_term_ids(c),
        is_lof_nagnag=c.lof_flags == NAGNAG_SITE,
        lof_filter_ids=_lof_filter_ids(c),
    )


def transcript_consequences_sort(
    ht: hl.Table,
) -> Callable[[hl.StructExpression], hl.StructExpression]:
    return lambda c: hl.bind(
        lambda is_coding, is_most_severe, is_canonical: (
            hl.cond(
                is_coding,
                hl.cond(
                    is_most_severe,
                    hl.cond(is_canonical, 1, 2),
                    hl.cond(is_canonical, 3, 4),
                ),
                hl.cond(
                    is_most_severe,
                    hl.cond(is_canonical, 5, 6),
                    hl.cond(is_canonical, 7, 8),
                ),
            )
        ),
        c.biotype_id == PROTEIN_CODING_ID,
        hl.set(c.consequence_term_ids).contains(
            TRANSCRIPT_CONSEQUENCE_TERMS_LOOKUP[ht.vep.most_severe_consequence],
        ),
        hl.or_else(c.canonical, 0) == 1,
    )
