---
name: Bad markdown
description: "incorrectly formatted outputs"
inputs:
    - ints1:
        - IntSequence1 | IntSequence2
        - list
parameters:
    - int1:
        - Int
        - int
outputs:
    - concatenated_ints:
        - IntSequence1
        - list
        - Oops
---
