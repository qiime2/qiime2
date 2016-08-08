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
