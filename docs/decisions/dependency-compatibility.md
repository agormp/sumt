# phylotreelib dependency compatibility

Date: 2026-10-22. Status: recommended policy; package requirements and version numbers have not been changed.

## Recommended policy

Use the oldest tested version that supplies all required features as the lower bound, and exclude the next incompatible major version. For example:

```text
phylotreelib>=2.3.1,<3
```

For sumt using the proposed phylotreelib 3.0 API, after that version exists:

```text
phylotreelib>=3.0.0,<4
```

The lower bound prevents installation with a library missing required features. The upper bound prevents automatically upgrading to a major release with incompatible API changes. A blanket `2.*` constraint alone would still permit versions older than required. `~=2.3.1` would restrict updates to the 2.3 series, not all compatible 2.x releases, so explicit bounds are clearer here.

This relies on maintaining the public API throughout a major series: compatible fixes use patch versions, compatible additions use minor versions, and incompatible public changes require a major version. Document which APIs are public, add deprecation warnings before removals where practical, and test sumt against both the minimum supported phylotreelib and the latest supported release. Bounds are a compatibility policy, not a substitute for tests.

Exact pins belong in a reproducible analysis environment or lock file, not necessarily sumt's general package requirement. The two projects' major numbers need not match: sumt can depend on phylotreelib 3 without becoming sumt 5 if its own supported user-facing interface remains compatible.

## Current development situation

The new `nexus(rooted=...)` argument is an additive library feature, but sumt calling it needs a newer minimum than published phylotreelib 2.3.1. The separate printing API removals are breaking changes planned for 3.0.0. An upper bound alone cannot supply a feature missing from an older library.

The source checkouts intentionally still carry their existing version numbers under Anders's instruction. Consequently, installed version metadata cannot distinguish these unreleased checkouts from the published packages. Test matched source checkouts together; do not publish them under the unchanged version numbers.

Before publishing: choose the library release version, update sumt's bounds to that version and the next major, test the actual built distributions in a fresh environment, publish the compatible library first, then publish sumt. Recheck compatibility if either project changes during that process. None of these release actions is currently authorized.

## Sources

Accessed 2026-10-22:

- [Python Packaging version-specifier specification](https://packaging.python.org/en/latest/specifications/version-specifiers/): semantics of comma-separated constraints, prefix matching, and compatible-release specifiers.
- [Semantic Versioning 2.0.0](https://semver.org/): compatible versus incompatible public API changes and major/minor/patch versioning.
