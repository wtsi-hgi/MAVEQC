# MAVEQC: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

Note that this changelog starts with version 0.3.0.229, as no log was maintained prior to this.

## 0.3.0.229 - [26th August 2026]

### Changed
* Library reads threshold for plasmid QC updated to 50% (screen QC library reads threshold unchanged)
* `sqc_library_per` in config.yaml replaced by separate screen (`screen_sqc_library_per`) and plasmid (`plasmid_sqc_library_per`) variables

### Fixed
* Fixed Docker-in-Docker settings to work with wider GitLab runner updates