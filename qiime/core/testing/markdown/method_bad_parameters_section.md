---
name: Bad markdown
description: "incorrectly formatted parameters"
inputs:
    - ints1:
        - IntSequence1 | IntSequence2
        - list
parameters:
    - int1:
        - Int
        - int
        - Oops
outputs:
    - concatenated_ints:
        - IntSequence1
        - list
---
