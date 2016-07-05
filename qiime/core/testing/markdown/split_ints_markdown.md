---
name: Split sequence of integers in half
description: "This method splits a sequence of integers in half, returning the two halves (left and right). If the input sequence's length is not evenly divisible by 2, the right half will have one more element than the left."
inputs:
    - ints:
        - IntSequence1
        - list
parameters: []
outputs:
    - left:
        - IntSequence1
        - list
    - right:
        - IntSequence1
        - list
---
### Find midpoint

```python
>>> middle = int(len(ints) / 2)
```

### Slice halves

```python
>>> left = ints[:middle]
>>> right = ints[middle:]
```
