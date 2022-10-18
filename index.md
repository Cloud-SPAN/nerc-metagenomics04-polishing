---
layout: lesson
root: .
---
In our previous lesson we used long reads to generate our draft assembly. Long read data is great for tasks like this because they produce a less fragmented assembly and are more likely to span areas with repeats. However, they are also more likely to contain sequencing errors than short read data.

We must therefore use further tools to improve the quality of our draft assembly. We can "polish" our assembly using both long and short read data. After that, we can perform quality control (QC) checks to see what impact the polishing has had.

By the end of this lesson you will be able to:
- polish a draft metagenome assembly using tools such as Medaka and Pilon
- check the quality of your draft assembly using Seqkit and metaQUAST