# phylotreelib dependency compatibility

Date: 2026-10-28. Updated: 2026-10-29 for phylotreelib 3.0.1 and sumt 4.2.1.

## Dependency policy

Use the oldest tested version that supplies all required features as the lower bound, and exclude the next incompatible major version. For sumt 4.2.1:

```text
phylotreelib>=3.0.1,<4
```

The lower bound prevents installation with a library missing required features. The upper bound prevents automatically upgrading to a major release with incompatible API changes. A blanket `3.*` constraint alone would still permit versions older than a later required minimum. `~=3.0.0` would restrict updates to the 3.0 series, not all compatible 3.x releases, so explicit bounds are clearer here.

This relies on maintaining the public API throughout a major series: compatible fixes use patch versions, compatible additions use minor versions, and incompatible public changes require a major version. Document which APIs are public, add deprecation warnings before removals where practical, and test sumt against both the minimum supported phylotreelib and the latest supported release. Bounds are a compatibility policy, not a substitute for tests.

Exact pins belong in a reproducible analysis environment or lock file, not necessarily sumt's general package requirement. The two projects' major numbers need not match: sumt 4.2.0 can depend on phylotreelib 3 while retaining its supported command-line interface.

## Release preparation

The patch release versions are phylotreelib 3.0.1 and sumt 4.2.1. The new minimum includes the Tree iteration and cache fixes. Both require Python 3.10 or newer.

phylotreelib 3.0.0 removes public printing settings and changes the input/output API. Rootedness is stored in `Tree.rooted`. Nexus output uses `print_rootedness` to include or suppress its declaration. The former `nexus(rooted=...)` argument has been removed. sumt 4.2.0 uses this API and cannot run against phylotreelib 2.3.1.

Before publishing, verify the actual built distributions together in fresh environments, including Python 3.10 and the current supported Python version. Check both the minimum supported phylotreelib and the latest supported release when they differ. Publish the compatible library first, then publish sumt. Recheck compatibility if either project changes during that process.

## Sources

Accessed 2026-10-22:

- [Python Packaging version-specifier specification](https://packaging.python.org/en/latest/specifications/version-specifiers/): semantics of comma-separated constraints, prefix matching, and compatible-release specifiers.
- [Semantic Versioning 2.0.0](https://semver.org/): compatible versus incompatible public API changes and major/minor/patch versioning.
