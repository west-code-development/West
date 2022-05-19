# Contributing to the source

Contributions are welcomed via merge requests. Contact the **WEST** developers before starting work to ensure it meshes well with the planned development direction and
standards set for the project.

## Version control

All changes in a pull request should be closely related. Multiple change sets that are loosely coupled should be proposed in separate pull requests. Use a consistent style for writing code.

## Features

New features should be applicable to a variety of use-cases. The **WEST** developers can assist you in designing flexible interfaces.

## Testing

Add tests for all new functionality.

## Release

We use [semantic versioning](https://semver.org/), i.e. version labels have the form v`<major>`.`<minor>`.`<patch>`

 - Patch release: v0.0.0 to v0.0.1, only bug fixes
 - Minor release: v0.0.0 to v0.1.0, bug fixes and new features that maintain backwards compatibility
 - Major release: v0.0.0 to v1.0.0, bug fixes and new features that break backwards compatibility

# Contributing to the documentation

Comment complex sections of code so that other developers can understand them.
Add demonstrations of new functionality, e.g. using Jupyter notebooks.

---

# Developer's workflow

  1. Clone the repository to your local machine
```bash
  $ git clone <http-git-repo>
```

  2. Use the following to see all the available branches and which branch you are on:
```bash
  $ git branch -a
```

  3. Create a branch for the feature you are going to work on. Internal users should branch from the `develop` branch.
```bash
  $ git checkout -b <myfeature>
```

  4. Make changes to the repository:
```bash
  $ touch <file>
```

  5. Check which files have been changed in your local directory by
```bash
  $ git status
```

  6. Commit your changes and push to the remote repository:
```bash
  $ git add <file>
  $ git commit -m "explain your modification"
  $ git push origin <myfeature>
```

  7. (Optional) Check commit history by
```bash
  $ git log
```

  8. Then create a merge request on the gitlab website that hosts the remote repository:

  9. You can do one of the following after the merge request has been approved:
     - Delete the branch `myfeature`
     - Keep the original branch and keep working on it

## Often-countered scenarios

  1. When you want to add another feature, create a new branch for the new feature, e.g. `another_feature`, branching from the latest `develop` branch.

  First switch to the `develop` branch:
```bash
  $ git checkout develop
```
  Then update the local `develop` branch to include changes to the remote:
```bash
  $ git pull origin develop
```
  Then create your new branch:
```bash
  $ git branch -b <another_feature>
```
  Then repeat #4-#9 above.

  2. Suppose you are working on your feature branch and other people have changed the `develop` branch in the remote repository after you created your `myfeature` branch, you need to merge the `develop` branch into your `myfeature` branch before submitting a merge request.

  To merge the `develop` branch into your `myfeature` branch, you need checkout your `myfeature` branch.
```bash
  $ git checkout <myfeature>
```
  Then do the following to merge the `develop` branch into the branch you are on it
```bash
  $ git pull origin develop
```
  Then do #6-#9 above.
