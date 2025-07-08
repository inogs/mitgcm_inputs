## MITgcm input generator

Welcome to the MITgcm input generator repository! This repository contains
tools that can be used to generate the static files that describe the
physical and biogeochemical properties of a geographical domain.
These files are specifically designed for data that has been discretized
using the [`bathytool`](https://github.com/inogs/bathytools) software.

**Note:**
This repository uses git LFS; please remember to execute
```bash
git lfs pull
```
after having cloned it.

### Instructions for Developers
This project uses **Poetry** as its package manager. To set up your
development environment, follow these steps:

1. Navigate to the project's root directory and execute the following
   command to install the required dependencies and create a virtual
   environment:
   ```bash
   poetry install
   ```

2. After setting up the environment, enable the `git` hooks by running:
   ```bash
   poetry run pre-commit install
   ```

These hooks ensure proper code formatting and linting during development.
