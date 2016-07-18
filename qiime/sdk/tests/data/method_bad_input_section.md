---
name: Bad markdown
description: "incorrectly formatted inputs"
inputs:
    - ints1:
        - IntSequence1 | IntSequence2
        - list
        - Oops
parameters:
    - int1:
        - Int
        - int
outputs:
    - concatenated_ints:
        - IntSequence1
        - list
---
