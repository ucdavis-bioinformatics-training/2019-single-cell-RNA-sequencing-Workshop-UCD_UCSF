# Running jobs on the cluster and using modules

Nik Joshi

najoshi@ucdavis.edu


**1\.** Here at the UC Davis Bioinformatics Core we have a large computational cluster (named lssc0) that we use for our analyses.

** We aren't going to use the cluster (just the workshop login node) for this workshp. **

**2\.** However, let's look at the 'module' command and its sub-commands. You will NOT find the 'module' command on all linux computers. Using modules is generally something that is used on a cluster and generally installed by a system administrator. The module system allows you to easily access any piece of software on the system from any directory. It basically changes your PATH variable (and possibly other environment variables) so that the shell searches the correct directories for the software you want to use. First, take a look at all the software available on our system:

    module avail

This is a list of all the software (with different versions) that you can access. Now try running the 'scythe' software:

    scythe

You should get an error saying that the command was not found. Take a look at your PATH variable.

    echo $PATH

These are the directories (colon separated) that are searched for anything you run on the command-line. In order to access a piece of software that is not in one of these default directories, we need to use the 'module load' command:

    module load scythe
    scythe --help

Use the 'which' command to find out where the 'scythe' command is actually located:

    which scythe

You'll see that scythe is located in a completely different place and yet you are able to access it. This is because the module command changes your PATH variable so that it has the correct directory. Take a look at your PATH again:

    echo $PATH

You'll see that the directory for scythe has been added to PATH.

---

**6\.** You can also load older versions of software. When you load a module without specifying a version, it will load the default (generally the latest) version. If you need an older version, you need to add the version number:

    module load bwa/0.6.2

---

**7\.** A few more module sub-commands that are useful:

    module purge

'module purge' will unload all of your modules. Which simply means that it will take out the directories for all of the modules from your PATH variable. Take a look at $PATH now:

    echo $PATH

The scythe and bwa directories are gone.

    module load scythe
    module rm scythe

'module rm' will unload the one module you specify.

    module load scythe
    module load sickle
    module load bwa
    module list

Finally, 'module list' will list all of your currently loaded modules in this terminal/session.
