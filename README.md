# NumberTheoreticTransform

Contains information related to my master's thesis that I finished in 1995.

## Description

When I attended graduate school at Marquette University, I needed to pick a topic for a master's thesis. I finally landed on high-precision arithmetic algorithms, specifically one called the number theoretic transform, or NTT. This was quite the challenge for me, because I didn't know anything about rings, fields, and other abstract algebra concepts. I also wasn't great at C. But, after a lot of research and hacking, I was able to demonstrate a working implementation that calculated pi to 5,000 digits of precision. Granted, in today's age, 5,000 digits of pi is a paltry amount. In fact, when I was writing my thesis, the pi record at that time was just over 4 billion digits. The current record (as of 2026.01.26) is 314 **trillion** digits. But my code **did** get 5,000 digits correct, and it was such a relief when it finally did the right thing!

I've always wanted to have my own digital version of the thesis. Unfortunately, I didn't keep a copy of the C code or the Word Perfect document I used to write the thesis. Fortunately, I did keep a hard copy of my thesis. I've always had the idea of somehow taking that content and getting it into a GitHub repository, but I also knew it would be a lot of tedious work to do that. But recently, Matt Parker has created the "Moon Pi" project, where an upcoming moon lander mission will be used in a creative way to calculate Pi on the moon! I've backed the project such that I could get my own pi estimate if all goes well, and...that's motivated me to finally get my own approach online.

## Goals

Here's a list of what I want to include in this repository:

* An approximation of the C code in C#. The C code isn't good. I know it isn't. But I don't want to support a version in C either. Therefore, my first version will be a faithful translation of the C code into C#. There will be no attempt to optimize or prettify it; I just want it to work.
* A modernized, improved version in C#. Once I get the original code in working order, I'll then try to update and optimize it as much as I can. I won't change the original, I'll just make another implementation based on the original and refactor that.
* A digitized version of my thesis. I'll probably just scan the thesis, take the images, and stitch them into a PDF.
* A digitized version of all my notes. I also kept **all** my notes, e-mails with my advisor, etc. and I want to scan them into a PDF as well.

## Additional Documentation

* [Marquette Reference](https://epublications.marquette.edu/theses/3962/)
* [Number Theoretic Transform](https://en.wikipedia.org/wiki/Discrete_Fourier_transform_over_a_ring#Number-theoretic_transform)
* [Matt Parker's Moon Pi](https://pi.space/)
* [Chronology of computation of pi](https://en.wikipedia.org/wiki/Chronology_of_computation_of_pi)
* [Chudnovsky algorithm](https://en.wikipedia.org/wiki/Chudnovsky_algorithm)

## Feedback

If you run into any issues, please add them [here](https://github.com/JasonBock/NumberTheoreticTransform/issues).