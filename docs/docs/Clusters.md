Running on clusters
----------------------------------
There is a file in the `templates` called `cluster.yaml`. This can be used to modify ressources needed for your data. I generally recommand moving the file to the root of the folder so that it doesn't get replaced by updates.

Bellow is an example of running on a cluster using the template file `cluster.yaml` on SLURM.

```
snakemake --cluster 'sbatch -n {cluster.n}  -t {cluster.time} --clusters=CLUSTERNAME --output={cluster.output}' --jobs N --cluster-config cluster.yaml --use-conda --local-cores C
```

* N: is the number of jobs you are allowed to run at the same time
* C: is the local-cores of the host machine. A few simple rules are gonna be run locally (not sent to nodes) because they are not that heavy (mostly plotting)
* CLUSTERNAME: the name of the cluster you want to use

Note: The default path for cluster logs in the cluster.yaml is `logs/cluster/`. If that folder doesn't exist, our cluster can't write and will crash without an error message.