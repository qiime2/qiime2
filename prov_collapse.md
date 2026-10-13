# Overall

We want a way to reduce the amount of space on disk/inodes used by provenance. There is some provenance we don't care about, particularly provenance from internal map/reduce steps. For instance, when you run a boots action with 1000 resamples, we don't need to keep around the provenance from all 1000 bootstrapped tables.

Ideally, we want this solution to not break forwards compatibility with prior versions of QIIME 2. We also want this solution to be as simple as possible for users and pipeline developers. In a perfect world, it would be invisible for them (especially end users).

The end result would be the provenance for pipeline output results would only show the provenance for the pipeline run and contain no provenance for any actions running internal to the pipeline.

## Pros:

1. Reduces space/inodes used on disk
2. Allows faster loading/parsing of prov

## Cons:

1. Could cause issues with recycling depending on approach
2. Could break forwards compatibility
3. Need to ensure we don't destroy provenance we actually want

# Approaches:

## 1 Redact internal prov after the pipeline has completed

After the pipeline has returned and all terminal results are collated, redact and delete all provenance pertaining to what happened inside the pipeline.

This could be in the form of a `ctx.drop_inner_prov` method called by the developer at the end of the pipeline, or a similar method called by us in the framework after the pipeline executes. The later would do this automatically for all pipelines. The former would require the dev to do it.

This approach can, to some extent, be applied to any of the other approaches. By this I mean, we could write the provenance as normal and then redact it post facto for only certain sections of the pipeline, or only for pipelines the developer said shouldn't keep their prov around. etc. etc.

### Pros:

1. Allows us to recycle as normal because the intermediate results will be there unless and until the pipeline succeeds (if --no-recycle is set then we could avoid storing the prov data at all)
2. Neither users nor devs need to be aware of it

### Cons:

1. Still does have to write to disk before deleting, so will take up space/inodes for a while

## 2 Some kind of context manager allowing redaction of specified portions of prov within the pipeline

Something like the following preventing tracking of provenance for the action in the manager:

```{python}
with ctx.no_prov():
    bar = ctx.get_action('foo', 'bar')
    bar()
```

Any actions executed within the context manager, or perhaps any actions gotten within the context manager, will not have their internal provenance recorded (or kept around post execution anyway).

### Pros:

1. Allows for precise control over what does and does not get recorded
2. Can be dynamic. Only ignore prov if certain conditions met e.g. only ignore prov of resampled tables if n >= 100

### Cons:

1. Requires the dev to know what they're doing
2. Need to be careful to put all necessary pieces in context manager and nothing extra

## 3 User controlled context manager/arg to drop internal prov

User control over when provenance isn't recorded. Something like the following.

On CLI:

```{bash}
qiime plugin pipeline --no-prov
```

In Python:

```{python}
from qiime2.sdk import no_prov

with no_prov():
    my_pipeline()
```

If we want to give the user A LOT of control, we could allow them to parameterize which actions in the pipeline don't record/drop prov. This would get very fiddly for the user, and for us. If we did go this far, it would probably be in our best interest to also have the generic full pipeline option.

### Pros:

1. Allows for user to control when prov is or isn't recorded

### Cons:

1. Requires the user to actually know what they're doing.

## 4 Drop Provenance after call to context method

Everything after ctx.no_prov is called gets dropped, or never recorded.

```{Python}
def my_pipeline(ctx):
    ctx.no_prov()

    ...
```

We could require this be called at the top of the pipeline, or we could allow it to be called wherever and drop prov on the rest of the pipeline post call.

### Pros:

1. Simple to use
2. Can be parameterized similar to #2

### Cons:

1. The pipeline developer needs to know to use this, but this is probably easier to use than #2
2. No granular control
