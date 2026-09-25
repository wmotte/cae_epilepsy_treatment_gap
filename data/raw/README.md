# data/raw

`wpp2024_lifetable_extract.tsv` is committed. It is a small extract of the UN
World Population Prospects 2024 abridged life table (both sexes, 2023) for
Nigeria, Ecuador and the UK, and the pipeline reads it on a fresh checkout.

The full workbook, `WPP2024_MORT_F07_1_ABRIDGED_LIFE_TABLE_BOTH_SEXES.xlsx`, is
not committed because at 131 MB it is over GitHub's file limit. You can download
it from <https://population.un.org/wpp/> ("Abridged life table, both sexes"). If
you place it in this folder, `R/01_prepare_data.R` rebuilds the extract from it.

Source: United Nations, Department of Economic and Social Affairs, Population
Division. World Population Prospects 2024. Licensed under CC BY 3.0 IGO.
