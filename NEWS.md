# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.6.8] - 2024.08.29

### Added

- Added `SparseMatrixCSC` and `AbstractMatrix` constructors.

### Fixed

- Fixed show method for SparseMatrixCSR.
- Fixes related to `TransposeFactorization` (see https://github.com/JuliaLang/julia/pull/46874).

## [0.6.6] - 2021.11.26

### Added

- Implemented `Base.setindex!`.

## [0.6.5] - 2021.10.20

### Fixed

- Return value of `LinearAlbegra.fillstored!`.

### Added

- Implemented `LinearAlbegra.rmul!`.

## [0.6.4] - 2021.10.20

### Added

- Implemented `LinearAlbegra.fillstored!`.

*Previous releases are not included in this Changelog*
