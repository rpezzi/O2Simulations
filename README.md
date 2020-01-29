# O2Simulations

Scripts and macros for MFT Simulations

# Getting started with O2 Simulations

## Table of Contents

1.  [Setup tools](#orgf5c8ee3)
2.  [Relevant skills](#org0544175)
3.  [References](#org54c5d9c)
4.  [Proposed activities (Chalenges)](#orgc25d1eb)


<a id="orgf5c8ee3"></a>

# Setup tools

-   O2
    -   <https://alice-doc.github.io/alice-analysis-tutorial/building/>
-   Emacs with Babel (to get organized - Optional?)
    -   <http://www.cachestocaches.com/2018/6/org-literate-programming/>
    -   <https://orgmode.org/worg/org-contrib/babel/>


<a id="org0544175"></a>

# Relevant skills

-   Root / macros / TTree
    -   Root Primer: <https://root.cern.ch/root/htmldoc/guides/primer/ROOTPrimer.html>
    -   Meet a TTree: <https://root.cern.ch/meet-ttree>
    -   Running MFT Simulations: <https://alice-talk.web.cern.ch/t/instructions-to-run-mft-simulations-wiki-post/470>
    -   Reconstruction workflow:
        -   ![img](https://alice-offline.web.cern.ch/sites/alice-offline.web.cern.ch/files/images/Reconstruction-Framework.gif)
        -   <https://alice-offline.web.cern.ch/Activities/Reconstruction/index.html>
-   version control with GIT
    -   GIT Pro book
        -   <https://git-scm.com/book/en/v2>
    -   ALICE Git tutorial
        -   <https://alisw.github.io/git-tutorial/>


<a id="org54c5d9c"></a>

# References

-   O2 Project DOxigen documentation: <https://aliceo2group.github.io/AliceO2/index.html>
-   ROOT Reference Documentation: <https://root.cern.ch/doc/master/>
-   C++ reference: <https://en.cppreference.com/w/>
-   ALICE Talk (O2 Discussion Forum)

<https://alice-talk.web.cern.ch/>


<a id="orgc25d1eb"></a>

# Proposed activities (Chalenges)

-   Find optimal QED Background Ymin for MFT simulations
-   Fix MFTMultiplicityEstimatorFromClusters.C
    -   <https://github.com/AliceO2Group/AliceO2/commit/43d3cb9d33f21f3c50aad9c57934aeb81ad7558d>
    -   <https://github.com/AliceO2Group/AliceO2/commit/f2713edc020f1ca1a43ee4fd410a049d22e04388#diff-9db8edc429ee4bbbc1f7869f3431d246>
    -   <https://github.com/AliceO2Group/AliceO2/commit/cfa31c0bb470798269947b8b6ff5a956e122226e>
-   Macro to estimate tracking efficiency as a function of damaged MFT Chips or MFT ladders


