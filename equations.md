Equations for models
================
Lucas Nell
11/10/2017

-   [Model A](#model-a)
    -   [Solving unknown parameters](#solving-unknown-parameters)
    -   [Differential equations](#differential-equations)
-   [Model B](#model-b)
    -   [Solving for unknown parameters](#solving-for-unknown-parameters)
    -   [Differential equations](#differential-equations-1)

Model A
=======

Solving unknown parameters
--------------------------

$$
s\_N = i\_N - a\_{NP} N\_{eq} P\_{eq} (1-P\_{eq}/k\_P) + (1-l\_D) m\_D D\_{eq} \\\\
s\_D = (1-l\_P) m\_P P\_{eq} + (1-l\_V) m\_V V\_{eq} + (1-l\_H) m\_H H\_{eq} + 
    (1-l\_R) m\_R R\_{eq} - a\_{DV} D\_{eq} V\_{eq} (1-V\_{eq}/k\_V) - m\_D D\_{eq} \\\\
s\_P = a\_{NP} N\_{eq} P\_{eq} (1-P\_{eq}/k\_P) -  a\_{PH} P\_{eq} H\_{eq} (1-H\_{eq}/k\_H) - m\_P P\_{eq} \\\\
s\_V = a\_{DV} D\_{eq} V\_{eq} (1-V\_{eq}/k\_V) - (a\_R V\_{eq} R\_{eq}) (1-R\_{eq}/k\_R) - m\_V V\_{eq} \\\\
s\_H = a\_{PH} P\_{eq} H\_{eq} (1-H\_{eq}/k\_H) - (a\_R H\_{eq} R\_{eq}) (1-R\_{eq}/k\_R) - m\_H H\_{eq} \\\\
s\_R = (a\_R V\_{eq} R\_{eq} + a\_R H\_{eq} R\_{eq}) (1-R\_{eq}/k\_R) - m\_R R\_{eq} \\\\
$$

Differential equations
----------------------

$$
\\begin{aligned}
N &= i\_N - a\_{NP} N P (1-P/k\_P) + (1-l\_D) m\_D D \\\\
D &= (1-l\_P) m\_P P + (1-l\_V) m\_V V + (1-l\_H) m\_H H + (1-l\_R) m\_R R + (1 - l\_M) m\_M M - a\_{DV} D V (1-V/k\_V) - m\_D D \\\\
P &= a\_{NP} N P (1-P/k\_P) - a\_{PH} P H (1-H/k\_H) - m\_P P \\\\
V &= a\_{DV} D V (1-V/k\_V) - (a\_R V R) (1-R/k\_R) - m\_V V \\\\
H &= a\_{PH} P H (1-H/k\_H) - (a\_R H R) (1-R/k\_R) - m\_H H \\\\
R &= (a\_R V R + a\_R H R + \\frac{a\_R M R}{2} ) (1-R/k\_R) - m\_R R \\\\
M &= i\_M - m\_M M - (\\frac{a\_R M R}{2}) (1-R/k\_R) \\\\
\\end{aligned}
$$

Model B
=======

Solving for unknown parameters
------------------------------

$$
\\begin{aligned}
s\_N &= i\_N - a\_{NP} N\_{eq} P\_{eq} (1-P\_{eq}/k\_P) + (1-l\_D) m\_D D\_{eq} \\\\
s\_D &= (1-l\_P) m\_P P\_{eq} + (1-l\_V) m\_V V\_{eq} + (1-l\_H) m\_H H\_{eq} + (1-l\_R) m\_R R\_{eq} - a\_{DV} D\_{eq} V\_{eq}/(1+a\_{DV} h\_D D\_{eq}) - m\_D D\_{eq} \\\\
s\_P &= a\_{NP} N\_{eq} P\_{eq} (1-P\_{eq}/k\_P) - a\_{PH} P\_{eq} H\_{eq}/(1+a\_{PH} h\_P P\_{eq}) - m\_P P\_{eq} \\\\
s\_V &= a\_{DV} D\_{eq} V\_{eq}/(1+a\_{DV} h\_D D\_{eq}) - (a\_R V\_{eq} R\_{eq})/(1+a\_R h\_R (V\_{eq}+H\_{eq})) - m\_V V\_{eq} \\\\
s\_H &= a\_{PH} P\_{eq} H\_{eq}/(1+a\_{PH} h\_P P\_{eq}) - (a\_R H\_{eq} R\_{eq})/(1+a\_R h\_R (V\_{eq}+H\_{eq})) - m\_H H\_{eq} \\\\
s\_R &= (a\_R V\_{eq} R\_{eq} + a\_R H\_{eq} R\_{eq})/(1+a\_R h\_R (V\_{eq}+H\_{eq})) - m\_R R\_{eq} \\\\
\\end{aligned}
$$

Differential equations
----------------------

$$
\\begin{aligned}
N &= i\_N - a\_{NP} N P (1-P/k\_P) + (1-l\_D) m\_D D \\\\
D &= (1-l\_P) m\_P P + (1-l\_V) m\_V V + (1-l\_H) m\_H H + (1-l\_R) m\_R R + (1-l\_M) m\_M M - a\_{DV} D V/(1+a\_{DV} h\_D D) - m\_D D \\\\
P &= a\_{NP} N P (1-P/k\_P) - a\_{PH} P H/(1+a\_{PH} h\_P P) - m\_P P \\\\
V &= a\_{DV} D V/(1+a\_{DV} h\_D D) - (a\_R V R)/(1+a\_R h\_R (V+H+M)) - m\_V V \\\\
H &= a\_{PH} P H/(1+a\_{PH} h\_P P) - (a\_R H R)/(1+a\_R h\_R (V+H+M)) - m\_H H \\\\
R &= (a\_R V R + a\_R H R + a\_R M R)/(1+a\_R hR (V+H+M)) - mR R \\\\
M &= i\_M - m\_M M - (a\_R M R)/(1+a\_R h\_R (V+H+M)) \\\\
\\end{aligned}
$$
