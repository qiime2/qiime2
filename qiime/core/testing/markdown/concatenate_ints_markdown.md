---
name: Concatenate integers
description: "This method concatenates integers into a single list in the order
they are provided."
inputs:
    - ints1:
        - IntSequence1 | IntSequence2
        - list
    - ints2:
        - IntSequence1
        - list
    - ints3:
        - IntSequence2
        - list
parameters:
    - int1:
        - Int
        - int
    - int2:
        - Int
        - int
outputs:
    - concatenated_ints:
        - IntSequence1
        - list
---
## Concatenate some integers

This method concatenates integers in the following steps.

### Join lists together

```python
>>> concatenated_ints = ints1 + ints2 + ints3
```

### Append integers

```python
>>> concatenated_ints.append(int1)
>>> concatenated_ints.append(int2)
```
