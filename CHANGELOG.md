# Version 2.0.5 (2016-10-08)

* REF: Refactored Archiver logic (#165)

* ENH: Result.extract returns path to extracted directory (#161)

* ENH: Pluralize format and type registration

* ENH: Add explicit format registration

* ENH: adds support for None as default value (#158)

* ENH: ports and improves redirected_stdio from qiime-studio (#156)

* ENH: Plugin util to directly invoke a transformer (#157)

* ENH: Update temp directory name to use prefix (#154)

* ENH: Port TestPluginBase to framework (#152)

* ENH/REF: add PluginManager.iter_entry_points(), add qiime.sdk.util (#150)

* MAINT: typo fix (#149)

* ENH: added parameter defaults to qiime.core.type.signature classes (#125)

* ENH: Archive version 0.3.0: use archive UUID as zipfile root directory (#142)

* ENH: ignore dotfiles in Archiver.save() (#137)

* ENH: Archiver.save() filepath validation (#136)

* ENH: Artifact/Visualization.save() appends file extension automatically (#135)

* ENH: Actions return Results object regardless of number of outputs (#131)

* ENH: Make formats smarter (#133)

* ENH: Artifact/Visualization __eq__/__ne__ (#130)

* ENH: Artifacts are coerced by transformers (#120)

* REF: Method/Visualizer -> Action (#118)

# Version 2.0.2 (2016-08-08)

* Added support for QIIME 1 feature metadata files to `qiime.Metadata`.

* Added `PluginManager.get_data_layout` for retrieving a `DataLayout` for a given semantic type.

* Improved formatting of .qza/.qzv archive metadata (`metadata.yml`).

* Added recommended archive format extensions to `Artifact` and `Visualization`.

* Added commonly referenced URLs to `qiime.sdk`.

* Added `qiime.plugin.plugin_init` to initialize a plugin package from a [cookiecutter](https://cookiecutter.readthedocs.io/en/latest/) template.

* Added `user_support_text` and `citation_text` to `qiime.plugin.Plugin`.

# Version 2.0.1 (2016-07-14)

Initial alpha release. At this stage, major backwards-incompatible API changes are expected to happen.
