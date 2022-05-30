- There are two ways to generate pdf file:

1. GitHub Action

   GitHub will automatically compile the paper each time the repository is updated. The pdf is available via the Actions tab in the project and click on the latest workflow run.

2. Docker
```
docker run --rm \
    --volume $PWD:/data \
    --user $(id -u):$(id -g) \
    --env JOURNAL=joss \
    openjournals/paperdraft
```

- [Official documentation](https://joss.readthedocs.io/en/latest/submitting.html#what-should-my-paper-contain) on what should be included in the paper.

