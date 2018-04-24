# Chapter 11: Fitting birth-death models

[pdf version]({{ site.baseurl }}/pdf/chapter11_fitbd.pdf)

## Section 11.1: Introduction, diversity hotspots

The number of species on the Earth remains highly uncertain, but our best estimates are around 10 million [@Mora2011-jh]. This is a mind-boggling number, and far more than have been described. As far as we know, all of those species are descended from a single common ancestor that lived some 4.2 billion years ago [@Hedges2009-lj]. All of these species formed by the process of speciation, the process by which one species splits into two (or more) descendants [@Coyne2004-vv].

Some parts of the tree of life have more species than others. This imbalance in diversity tells us that speciation is much more common in some lineages than others [@Mooers1997-ow]. Likewise, numerous studies have argued that certain habitats are "hotbeds" of speciation [e.g. @Hutter2017-tf, @Miller2017-ja]. For example, the high Andes ecosystem called the Páromo - a peculiar landscape of alien-looking plants and spectacled bears - might harbor the highest speciation rates on the planet [@Madrinan2013-oi].


![Figure 11.1 Páromo ecosystem, Chingaza Natural National Park, Colombia. Photo taken by the author.]({{ site.baseurl }}/images/figure11-1.png)

In this chapter we will explore how we can learn about speciation and extinction rates from the tree of life. We will use birth-death models, simple models of how species form and go extinct through time. Birth-death models can be applied to data on clade ages and diversities, or fit to the branching times in phylogenetic trees. We will explore both maximum likelihood and Bayesian methods to do both of these things.

Key questions:

1. How can we calculate speciation and extinction rates using clade ages and diversities?
2. How can we fit a birth-death model to the pattern of branching times on a phylogenetic tree?

## Section 11.2: Clade age and diversity

If we know the age of a clade and its current diversity, then we can calculate the net diversification rate for that clade. Before presenting equations, I want to make a distinction between two different ways that one can measure the age of a clade: the stem age and the crown age. A clade's crown age is the age of the common ancestor of all species in the clade. By contrast, a clade's stem age measures the time that that clade descended from a common ancestor with its sister clade. The stem age of a clade is always at least as old, and usually older, than the crown age.

![Figure 11.2. Crown and stem age of a clade.]({{ site.baseurl }}/images/figure11-2.png)


Magallón and Sanderson [-@Magallon2001-xi] give an equation for the estimate of net diversification rate given clade age and diversity:

(eq. 11.1)

<div>
$$
\hat{r} = \frac{ln(n)}{t_{stem}}
$$
</div>

where $t_{stem}$ is the stem group age. In this equation we also see the net diversification rate parameter, $r=\lambda-\mu$ (see Chapter 10). This is a good reminder that the parameter that best predicts how species accumulate through time reflects the balance between speciation and extinction rates.

Alternatively, one can use $t_{crown}$, the crown group age:

(eq. 11.2)

<div>
$$
\hat{r} = \frac{ln(n)-ln(2)}{t_{crown}}
$$
</div>

The two equations differ because at the crown group age one is considering the clade's diversification starting with two lineages rather than one (Figure 11.2).

Even though these two equations reflect the balance of births and deaths through time, they give ML estimates for $r$ only under a pure birth model where there is no extinction. If there has been extinction in the history of the clade, then our estimates using equations 11.1 and 11.2 will be biased. The bias comes from the fact that we see only clades that survive to the present day, and miss any clades of the same age that died out before they could be observed. By observing only the "winners" of the diversification lottery, we overestimate the net diversification rate. If we know the relative amount of extinction, then we can correct for this bias.

Under a scenario with extinction, one can define $\epsilon= \mu / \lambda$ and use the following method-of-moments estimators from Rohatgi [-@Rohatgi1976-du, following the notation of @Magallon2001-xi]:

(eq. 11.3)

<div>
$$
\hat{r} = \frac{log⁡[n(1-\epsilon)+\epsilon]}{t_{stem}}
$$
</div>

for stem age, and

(eq. 11.4)

<div>
$$
\hat{r} = \frac{ln[\frac{n(1-\epsilon^2 )}{2}+2 \epsilon+\frac{(1-\epsilon) \sqrt{n(n \epsilon^2-8\epsilon+2 n \epsilon+n)}}{2}]-ln(2)}{t_{crown}}
$$
</div>

for crown age. (Note that eq. 11.3 and 11.4 reduce to 11.1 and 11.2, respectively, when $\epsilon = 0$). Of course, we usually have little idea what $\epsilon$ should be. Common practice in the literature is to try a few different values for $\epsilon$ and see how the results change [e.g. @Magallon2001-xi].

Magallón and Sanderson [-@Magallon2001-xi], following Strathmann and Slatkin [-@Strathmann1983-lz], also describe how to use eq. 10.13 and 10.15 to calculate confidence intervals for the number of species at a given time.

As a worked example, let's consider the data in table 11.1, which give the crown ages and diversities of a number of plant lineages in the Páromo [from @Madrinan2013-oi]. For each lineage, I have calculated the pure-birth estimate of speciation rate (from equation 11.2, since these are crown ages), and net diversification rates under three scenarios for extinction ($\epsilon = 0.1$,   $\epsilon = 0.5$, and $\epsilon = 0.9$).

| Lineage |	n |	Age	| $\hat{r}_{pb}$ | $\hat{r}_{\epsilon = 0.1}$ | $\hat{r}_{\epsilon = 0.5}$ | $\hat{r}_{\epsilon = 0.9}$ |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Aragoa | 17 | 0.42 | 5.10 | 5.08 | 4.53 | 2.15 |
| Arcytophyllum | 14 | 10.96 | 0.18 | 0.18 | 0.16 | 0.07 |
| Berberis | 32 | 3.80 | 0.73 | 0.73 | 0.66 | 0.36 |
| Calceolaria | 65 | 2.50 | 1.39 | 1.39 | 1.28 | 0.78 |
| Draba | 55 | 3.05 | 1.09 | 1.08 | 1.00 | 0.59 |
| Espeletiinae | 120 | 4.04 | 1.01 | 1.01 | 0.94 | 0.62 |
| Festuca | 36 | 4.28 | 0.68 | 0.67 | 0.61 | 0.34 |
| Jamesonia + Eriosorus | 32 | 7.60 | 0.36 | 0.36 | 0.33 | 0.18 |
| Lupinus | 66 | 1.47 | 2.38 | 2.37 | 2.19 | 1.34 |
| Lysipomia | 27 | 8.96 | 0.29 | 0.29 | 0.26 | 0.14 |
| Oreobolus | 5 | 3.01 | 0.30 | 0.30 | 0.26 | 0.09 |
| Puya | 46 | 0.80 | 3.92 | 3.91 | 3.58 | 2.07 |
| Valeriana | 53 | 14.58 | 0.22 | 0.22 | 0.21 | 0.12 |

Table 11.1. Estimates of net diversification rates for Páramo lineages [data from @Madrinan2013-oi] using equations 11.2 and 11.4.

Inspecting the last four columns of Table 11.1, we can make a few general observations. First, these plant lineages really do have remarkably high diversification rates [@Madrinan2013-oi]. Second, the net diversification rate we estimate depends on what we assume about relative extinction rates ($\epsilon$). You can see in the table that the effect of extinction is relatively mild until extinction rates are assumed to be quite high ($\epsilon = 0.9$). Finally, assuming different levels of extinction affects the diversification rates but not their relative ordering. In all cases, net diversification rates for Aragoa, which formed 17 species in less than half a million years, is higher than the rest of the clades. This relationship holds only when we assume relative extinction rates are constant across clades, though. For example, the net diversification rate we calculate for Calceolaria with $\epsilon = 0$ is higher than the calculated rate for Aragoa with $\epsilon = 0.9$. In other words, we can't completely ignore the role of extinction in altering our view of present-day diversity patterns.

We can also estimate birth and death rates for clade ages and diversities using ML or Bayesian approaches. We already know the full probability distribution for birth-death models starting from any standing diversity $N(0)=n_0$ (see equations 10.13 and 10.15). We can use these equations to calculate the likelihood of any particular combination of $N$ and $t$ (either $t_{stem}$ or $t_{crown}$) given particular values of $\lambda$ and $\mu$. We can then find parameter values that maximize that likelihood. Of course, with data from only a single clade, we cannot estimate parameters reliably; in fact, we are trying to estimate two parameters from a single data point, which is a futile endeavor. (It is common, in this case, to assume some level of extinction and calculate net diversification rates based on that, as we did in Table 11.1 above).

One can also assume that a set of clades have the same speciation and extinction rates and fit them simultaneously, estimating ML parameter values. This is the approach taken by Magallón and Sanderson's [-@Magallon2001-xi] paper on diversification rates across angiosperms. When we apply this approach to the Paramo data, shown above, we obtain ML estimates of $\hat{r} = 0.27$ and $\hat{\epsilon} = 0$. If we were forced to estimate an overall average rate of speciation for all of these clades, this might be a reasonable estimate. However, the table above also suggests that some of these clades are diversifiying faster than others. We will return to the issue of variation in diversification rates across clades in the next chapter.

We can also use a Bayesian approach to calculate posterior distributions for birth and death rates based on clade ages and diversities. This approach has not, to my knowledge, been implemented in any software package, although the method is straightforward [for a related approach, see @Hohna2016-dw]. To do this, we will modify the basic algorithm for Bayesian MCMC (see Chapter 2) as follows:

1.  Sample a set of starting parameter values, $r$ and $\epsilon$, from their prior distributions. For this example, we can set our prior distribution for both parameters as exponential with a mean and variance of $\lambda_{prior}$ (note that your choice for this parameter should depend on the units you are using, especially for $r$). We then select starting r and $\epsilon$ from their priors.

2.  Given the current parameter values, select new proposed parameter values using the proposal density $Q(p'|p)$. For both parameter values, we can use a uniform proposal density with width $w_p$, so that $Q(p'|p) ~ U(p-w_p/2,p+w_p/2)$. We can either choose both parameter values simultaneously, or one at a time (the latter is typically more effective).

3. Calculate three ratios:

a. The prior odds ratio. This is the ratio of the probability of drawing the parameter values $p$ and $p'$ from the prior. Since we have exponential priors for both parameters, we can calculate this ratio as:

(eq. 11.5)

<div>
$$
R_{prior} = \frac{\lambda_{prior} e^{-\lambda_{prior} p'}}{\lambda_{prior} e^{-\lambda_{prior} p}}=e^{\lambda_{prior} (p-p')}
$$
</div>

b. The proposal density ratio. This is the ratio of probability of proposals going from $p$ to $p'$ and the reverse. We have already declared a symmetrical proposal density, so that $Q(p'|p) = Q(p|p')$ and $R_{proposal} = 1$.

c. The likelihood ratio. This is the ratio of probabilities of the data given the two different parameter values. We can calculate these probabilities from equations 10.13 or 10.16 (depending on if the data are stem ages or crown ages).

4. Find $R_{accept}$ as the product of the prior odds, proposal density ratio, and the likelihood ratio. In this case, the proposal density ratio is 1, so:

(eq. 11.6)

<div>
$$
R_{accept} = R_{prior} \cdot R_{likelihood}
$$
</div>



5.	Draw a random number $u$ from a uniform distribution between 0 and 1. If $u < R_{accept}$, accept the proposed value of both parameters; otherwise reject, and retain the current value of the two parameters.

6.	Repeat steps 2-5 a large number of times.

When we apply this technique to the Páromo [from @Madrinan2013-oi] and priors with $\lambda_{prior} = 1$ for both $r$ and $\epsilon$, we obtain posterior distributions: $r$ (mean = 0.497, 95% CI = 0.08-1.77) and $\epsilon$ (mean = 0.36, 95% CI = 0.02-0.84; Figure 11.3). Note that these results are substantially different than the ML estimates. This is because in our Bayesian analysis our prior on extinction rates is exponential, giving a high probability relatively strong levels of extinction relative to speciation.

![Figure 11.3. Posterior distribution for $r$ and $\epsilon$ for Páromo clades [@Madrinan2013-oi].]({{ site.baseurl }}/images/figure11-3.png)

Thus, we can estimate diversification rates from data on clade ages and diversities. If we have a whole set of such clades, we can (in principal) estimate both speciation and extinction rates, so long as we are willing to assume that all of the clades share equal diversification rates. However, as we will see in the next section, this assumption is almost always dubious!

##  Section 11.3: Tree Balance

As we discussed in Chapter 10, tree balance considers how "balanced" the branches of a phylogenetic tree are. That is, if we look at each node in the tree, are the two sister clades of the same size (balanced) or wildly different (imbalanced)?

Birth-death trees have a certain amount of "balance," perhaps a bit less than your intuition might suggest (see chapter 10). We can look to real trees to see if the amount of balance matches what we expect under birth-death models. A less balanced pattern in real trees would suggest that speciation and/or extinction rate vary among lineages more than we would expect. By contrast, more balanced trees would suggest more even and predictable diversification across the tree of life than expected under birth-death models. This approach traces back to Raup and colleagues, who applied stochastic birth-death models to paleontology in a series of influential papers in the 1970s [e.g. @Raup1973-xy, @Raup1974-dx].

##  Section 11.3a: Sister clades and the balance of individual nodes

For single nodes, we already know that the distribution of sister taxa species richness is uniform over all possible divisions of $N_n$ species into two clades of size $N_a$ and $N_b$ (Chapter 11). This idea idea leads to simple test of whether the distribution of species between two sister clades is unusual compared to the expectation under a birth-death model [@Slowinski1993-ks]. This test can be used, for example, to test whether the diversity of exceptional clades, like passerine birds, is higher than one would expect when compared to their sister clade. This is the simplest measure of tree balance, as it only considers one node in the tree at a time.

 Slowinsky and Guyer [-@Slowinski1993-ks] developed a test based on calculating a P-value for a division at least as extreme as seen in a particular comparisons of sister clades. We consider $N_n$ total species divided into two sister clades of sizes $N_a$ and $N_b$, where $N_a < N_b$ and $N_a + N_b = N_n$. Then:

(eq. 11.7)

If $N_a \neq N_b$:
<div>
$$
P = \frac{2 N_a}{N_n - 1}
$$
</div>

If $N_a = N_b$ or $P > 1$ then set $P = 1$

For example, we can assess diversification in the Andean representatives of the legume genus *Lupinus* [@Hughes2006-nr]. This genus includes one young radiation of 81 Andean species, spanning a wide range of growth forms. The likely sister clade to this spectacular Andean radiation is a clade of *Lupinus* species in Mexico that includes 46 species [@Drummond2012-zs]. In this case $N_a = 81 - 46 = 35$, and we can then calculate a P-value testing the null hypothesis that both of these clades have the same diversification rate:

(eq. 11.8)

<div>
$$
P = \frac{2 N_a}{N_n - 1} = \frac{2 \cdot 35}{81 - 1} = 0.875
$$
</div>

We cannot reject the null hypothesis. Indeed, later work suggests that the actual increase in diversification rate for Lupinus occurred deeper in the phylogenetic tree, in the ancestor of a more broadly ranging New World clade [@Drummond2012-zs, @Hughes2006-nr].

Often, we are interested in testing whether a particular trait - say, dispersal into the Páramo - is responsible for the increase in species richness that we see in some clades. In that case, a single comparison of sister clades may be unsatisfying, as sister clades almost always differ in many characters, beyond just the trait of interest. Even if the clade with our putative "key innovation" is more diverse, we still might not be confidence in inferring a correlation from a single observation.

To address this problem, many studies have used natural replicates across the tree of life, comparing the species richnesses of many pairs of sister clades that differ in a given trait of interest. Following Slowinsky and Guyer's [-@Slowinski1993-ks] logic above, we could calculate a p-value for each clade, and then combine those p-values into an overall test. In this case, one clade (with diversity $N_1$) has the trait of interest and the other does not ($N_0$), and our formula is half of equation 11.5 since we will consider this a one-tailed test:

(eq. 11.9)

<div>
$$
P = \frac{N_0}{N_n - 1}
$$
</div>

When analyzing replicate clade comparisons - e.g. many sister clades, where in each case one has the trait of interest and the other does not - Slowinsky and Guyer [-@Slowinski1993-ks] recommended combining these p-values using Fisher's combined probability test, so that:

(eq. 11.10)

<div>
$$
{\chi^2}_{combined} = -2 \sum ln(P_i)
$$
</div>

Here, the $P_i$ values are from $i$ independent sister clade comparisons, each using equation 11.9. Under the null hypothesis where the character of interest does not increase diversification rates, the test statistic, ${\chi^2}_{combined}$, should follow a chi-squared distribution with $2k$ degrees of freedom where k is the number of tests. But before you use this combined probability approach, see what happens when we apply it to a real example!

As an example, consider the following data, which compares the diversity of many sister pairs of plants. In each case, one clade has fleshy fruits and the other dry [data from @Vamosi2005-mn]:

| Fleshy fruit clade | $n_{fleshy}$ |	Dry fruit clade	| $n_{dry}$ |
| --- | ---: | --- | ---: |
| *Pangium* | 1 | *Acharia*+*Kigellaria* | 2 |
| *Cyrilla* | 1 | *Clethra* | 64 |
| *Roussea* | 1 | *Lobelia* | 300 |
| *Myriophylum+Haloragis+Penthorum* | 1 | *Tetracarpaea* | 89 |
| *Austrobaileya* | 1 | *Illicium+Schisandra* | 67 |
| *Davidsonia* | 3 | *Bauera* | 4 |
| *Mitchella* | 3 | *Pentas* | 34 |
| *Milligania * | 5 | *Borya* | 10 |
| *Sambucus* | 9 |  *Viburnum* | 150 |
| *Pereskia* | 16 | *Mollugo* | 35 |
| *Decaisnea+Sargentodoxa+ Tinospora+Menispermum +  Nandina+Caulophyllum+Hydrastis+Glaucidium* | 33 | *Euptelea* | 2 |
| *Tetracera* | 40 | *Dillenia* | 60 |
| *Osbeckia* | 50 | *Mouriri* | 81 |
| *Hippocratea* | 100 | *Plagiopteron* | 1 |
| *Cyclanthus+Sphaeradenia+Freycinetia* | 216 | *Petrosavia+Japonlirion* | 3 |
| *Bixa* | 393 | *Theobroma+Grewia+Tilia+ Sterculia+Durio* | 1 |
| *Impatiens* | 850 | *Idria* | 11 |
| *Lamium+Clerodendrum+Callicarpa+Phyla +Pedicularis+Paulownia* | 947 | *Euthystachys* | 1 |
| *Callicarpa+Phyla+Pedicularis+Paulownia+Solanum* | 1700 | *Solanum* | 18 |

The individual clades show mixed support for the hypothesis, with only 7 of the 18 comparisons showing higher diversity in the fleshy clade, but 6 of those 7 comparisons significant at $P < 0.05$ using equation 11.9. The combined probability test gives a test statistic of ${\chi^2}_{combined} = 72.8$. Comparing this to a $\chi^2$ distribution with 36 degrees of freedom, we obtain $P = 0.00027$, a highly significant result. This implies that fleshy fruits do, in fact, result in a higher diversification rate.

However, if we test the opposite hypothesis, we see a problem with the combined probability test of equation 11.10. First, notice that 11 of 18 comparisons show higher diversity in the non-fleshy clade, with 4 significant at $P < 0.05$. The combined probability test gives ${\chi^2}_{combined} = 58.9$ and $P = 0.0094$. So we reject the null hypothesis and conclude that non-fleshy fruits diversify at a higher rate! In other words, we can reject the null hypothesis in both directions with this example.

What's going on here? It turns out that this test is very sensitive to outliers - that is, clades with extreme differences in diversity. These clades are very different than what one would expect under the null hypothesis, leading to rejection of the null - and, in some cases with two characters, when there are outliers on both sides [e.g. the proportion of species in each state has a u-shaped distribution; @Paradis2012-vp] we can show that both characters significantly increase diversity [@Vamosi2005-mn]!

Fortunately, there are a number of improved methods that can be used that are similar in spirit to the original Slowinsky and Guyer test but more statistically robust [e.g @Paradis2012-vp]. For example, we can apply the "richness Yule test" as described in Paradis [-@Paradis2012-vp], to the data from Vamosi et al. [-@Vamosi2005-mn]. This is a modified version of the McConway-Sims test [@McConway2004-jg], and compares the likelihood of a equal rate yule model applied to all clades to a model where one trait is associated with higher or lower diversification rates. This test requires knowledge of clade ages, which I don't have for these data, but [ @Paradis2012-vp] shows that the test is robust to this assumption and recommends substituting a large and equal age for each clade. I chose 1000, and found a significant likelihood ratio test (null model $lnL = -215.6$, alternative model $lnL = -205.7$, $P = 0.000008$). This method estimates a higher rate of diversification for fleshy fruits (since the age of the clade is arbitrary, the actual rates are not meaningful, but their estimated ratio $\lambda_1 / \lambda_0 = 1.39$ suggests that fleshy fruited lineages have a diversification rate almost 40% higher).

##  Section 11.3b: Balance of whole phylogenetic tremendous

We can assess the overall balance of an entire phylogenetic tree using tree balance statistics. As discussed, I will describe just one common statistic, Colless' I, since other metrics capture the same pattern in slightly different ways.

To calculate Colless' I, we can use equation 10.18. This result will depend strongly on tree size, and so is not comparable across trees of different sizes; to allow comparisons, $I_c$ is usually standardized by subtracting the expected mean for trees of that size under an ERM model, and dividing by the standard deviation. Both of these can be calculated analytically [@Blum2006-xo], and standardized I_c calculated using a small approximation [following @Bortolussi2006-jx] as:

(eq. 11.11)

<div>
$$
I^{'}_c = \frac{I_c-n*log(n)-n(\gamma-1-log(2))}{n}
$$
</div>

Since the test statistics are based on descriptions of patterns in trees rather than particular processes, the relationship between imbalance and evolutionary processes can be difficult to untangle! But all tree balance indices allow one to reject the null hypothesis that the tree was generated under a birth-death model. Actually, the expected patterns of tree balance are absolutely identical under a broader class of models called "Equal-Rates Markov" (ERM) models [@Harding1971-bb, @Mooers1997-ow]. ERM models specify that diversification rates (both speciation and extinction) are equal across all lineages for any particular point in time. However, those rates may or may not change through time. If they don't change through time, then we have a constant rate birth-death model, as described above - so birth-death models are ERM models. But ERM models also include, for example, models where birth rates slow through time, or extinction rates increase through time, and so on. All of these models predict exactly the same pattern of tree balance.

Typical steps for using tree balance indices to test the null hypothesis that the tree was generated under an ERM model are as follows:

1. Calculate tree balance using a tree balance statistic.
2. Simulate pure birth trees to general a null distribution of the test statistic. We are considering the set of ERM models as our null, but since pure-birth is simple and still ERM we can use it to get the correct null distribution.
3. Compare the actual test statistic to the null distribution. If the actual test statistic is in the tails of the null distribution, then your data deviates from an ERM model.

Step 2 is unnecessary in cases where we know null distributions for tree balance statistics analytically [e.g. @Blum2006-jo]. There are also some examples in the literature of considering null distributions other than ERM. For example, Mooers and Heard [-@Mooers1997-ow] consider two other null models, PDA and EPT, which consider different statistical distributions of tree shapes (but both of these are difficult to tie to any particular evolutionary process).

Typically, phylogenetic trees are more imbalanced than expected under the ERM model.  In fact, this is one of the most robust generalizations that one can make about macroevolutionary patterns in phylogenetic trees. This deviation means that diversification rates vary among lineages in the tree of life. We will discuss how to quantify and describe this variation in later chapters. These tests are all similar in that they use multiple non-nested comparisons of species richness in sister clades to calculate a test statistic, which is then compared to a null distribution, usually based on a constant-rates birth-death process [reviewed in @Vamosi2005-mn, @Paradis2012-vp].

As an example, we can apply the whole-tree balance approach to the tree of *Lupinus* [@Drummond2012-zs]. For this tree, which has 137 tips, we calculate $I_c = 1010$ and $I^{'}_c = 3.57$. This is much higher than expected by chance under an ERM model, with $P = 0.0004$. That is, our tree is significantly more imbalanced than expected under a ERM model, which includes both pure birth and birth-death. We can safely conclude that there is variation in speciation and/or extinction rates across lineages in the tree.

##  Section 11.4: Fitting birth-death models to branching times

Another approach that uses more of the information in a phylogenetic tree involves fitting birth-death models to the distribution of branching times. This approach traces all the way back to Yule [-@Yule1925-dv], who first applied stochastic process models to the growth of phylogenetic trees. More recently, a series of papers by Raup and colleagues [@Gould1977-dh, @Raup1974-dx, @Raup1985-wn, @Raup1973-xy] spurred modern approaches to quantitative macroevolution by simulating random clades, then demonstrating how variable such clades grown under simple birth-death models can be.

Most modern approaches to fitting birth-death models to phylogenetic trees use the intervals between speciation events on a tree - the "waiting times" between successive speciation - to estimate the parameters of birth-death models. Figure 10.2 shows these waiting times. Frequently, information about the pattern of species accumulation in a phylogenetic tree is summarized by a lineage-through-time (LTT) plot, which is a plot of the number of lineages in a tree against time (see Figure 10.9). As I introduced in Chapter 10, the y-axis of LTT plots is log-transformed, so that the expected pattern under a constant-rate pure-birth model is a straight line. Note also that LTT plots ignore the relative order of speciation events. Stadler [-@Stadler2013-vw] calls models justifying such an approach "species-exchangable" models - we can change the identity of species at any time point without changing the expected behavior of the model. Because of this, approaches to understanding birth-death models based on branching times are different from - and complementary to - approaches based on tree topology, like tree balance.

As discussed in the previous chapter, even though we often have no information about extinct species in a clade, we can still (in theory) infer the presence of extinction from an LTT plot. The signal of extinction is an excess of young lineages, which is seen as the "pull of the recent" in our LTT plots (Figure 10.10). Statistical approaches can capture this pattern in a more rigorous way.

###  Section 11.4a: Likelihood of waiting times under a birth-death model

In order to use ML and Bayesian methods for estimating the parameters of birth-death models from comparative data, we need to write down the likelihoods of the waiting times between speciation events in a tree. There is a little bit of variation in notation in the literature, so I will follow Stadler [@Stadler2013-vw] to maintain consistency. We will assume that the clade begins at time $t_1$ with a pair of species. Most analyses follow this convention, and condition the process as starting at the time $t_1$, representing the node at the root of the tree. This makes sense because we rarely have information on the stem age of our clade. We will also condition on both of these initial lineages surviving to the present day, as this is a requirement to obtain a tree with this crown age [e.g. @Stadler2013-yf equation 5).

Speciation and extinction events occur at various times, and the process ends at time $0$ when the clade has $n$ extant species - that is, we measure time backwards from the present day. Extinction will result in species that do not extend all the way to time 0. For now, we will assume that we only have data on extant species. We will refer to the phylogenetic tree that shows branching times leading to the extant species as the reconstructed tree [@Nee1994-xg]. For a reconstructed tree with $n$ species, there are $n-1$ speciation times, which we will denote as $t_1$, $t_2$, $t_3$, ..., $t_{n-1}$. The leaves of our ultrametric tree all terminate at time 0.

Note that in this notation, $t_1 > t_2 > \dots > t_{n-1} > 0$, that is, our speciation times are measured backwards from the tips, and as we increase the index the times are constantly decreasing (this is an important notational difference between Stadler [-@Stadler2013-yf], used here, and Nee [-@Nee1994-xg and others], the latter of which considers the time intervals between speciation events, e.g. $t_1 - t_2$ in our notation). For now, we will assume complete sampling; that is, all $n$ species alive at the present day are represented in the tree.

We will now derive a likelihood of of observing the set of speciation times $t_1$, $t_2$, ..., $t_{n-1}$ given the extant diversity of the clade, $n$, and our birth-death model parameters $\lambda$ and $\mu$. To do this, we will use an approach based on differential equations introduced by [xxx Maddison] and applied to the present context by [xxx FitzJohn et al.].

We need to keep track of two probabilities: $D_N(t)$, the probability that a lineage at some time t in the past will evolve into the extant clade N as observed today; and $E(t)$, the probability that a lineage at some time t will go completely extinct and leave no descendants at the present day. (Later, we will redefine $E(t)$ so that it includes the possibility that the lineages has descendants but none have been sampled). The general idea of how this method works is that we will assign values to these probabilities at the tips of the tree, and then define a set of rules to update them as we flow back through the tree from the tips to the root. When we arrive at the root of the tree, our $D_N(t)$ will represent the probability of observing the actual tree given our model - that is, the value of $D_N(t)$ at the root gives the likelihood.

So, what rules can we use to make these calculations backwards in time through the tree? First, we can define our starting points. Since every tip $i$ represents a living lineage, we can define $D_N(t) = 1$ and $E(t) = 0$ at each tip.

We now need to define how these probabilities change as we move backwards along branches, and what happens at the tree nodes.

Next, imagine we move backwards along some section of a tree branch with no nodes. Since that section of branch exists in our tree, we know two things: the lineage did not go extinct during that time, and if speciation occured, the lineage that split off did not survive to the present day. We can capture these two possibilities in a differential equation that considers how our overall liklihood changes over some very small unit of time.

(eq. 11.12)

<div>
$$
\frac{dD_N(t)}{dt} = -(\lambda + \mu) D_N(t) + 2 \lambda E(t) D_N(t)
$$
</div>

Here, the first part of the equation, $-(\lambda + \mu) D_N(t)$, represents the probability of not speciating nor going extinct, while the second part, $2 \lambda E(t) D_N(t)$,  represents the probability of speciation followed by the ultimate extinction of one of the two daughter lineages. The $2$ in this equation appears because we must account for the fact that, following speciation from an ancestor to daughters A and B, we would see the same pattern no matter which of the two descendants survived to the present.

We also need to update our extinction probability going back through the tree:

(eq. 11.13)

<div>
$$
\frac{dE(t)}{dt} = \mu - (\mu + \lambda) E(t) + \lambda E(t)^2
$$
</div>

The three parts of this equation represent the three ways a lineage might not make it to the present day: either it goes extinct during the interval being considered ($\mu$), it survives that interval but goes extinct some time later ($- (\mu + \lambda) E(t)$), or it speciates in the interval but both descendants go extinct before the present day ($\lambda E(t)^2$). We will also specify that $\lambda > \mu$; it is possible to relax that assumption, but it makes the solution more complicated.

We can solve these equations so that we will be able to update the probability moving backwards along any branch of the tree with length t. First, solving equation 11.13 and using our initial condition E(0) = 0:

(eq. 11.14)
<div>
$$
E(t) = 1 - \frac{\lambda-\mu}{\lambda - (\lambda-\mu)e^{(\lambda - \mu)t}}
$$
</div>


We can now substitute this expression for E(t) into eq. 11.12 and solve, conditioning on $D_N(0) = 1$:

(eq. 11.15)
<div>
$$
D_N(t) = e^{-(\lambda - \mu)(t - t_N)} \frac{(\lambda - (\lambda-\mu)e^{(\lambda - \mu)t_N})^2}{(\lambda - (\lambda-\mu)e^{(\lambda - \mu)t})^2} \cdot D_N(t_N)
$$
</div>

Here t_N is the depth (measured from the present day) of node N.

Now we need to consider what happens when two branches come together at a node. Since there is a node, we know there has been a speciation event. We multiply the probability calculations flowing down each branch by the probability of a speciation event. So:

(eq. 11.16)
<div>
$$
D_{N'}(t) = D_{N}(t) D_{M}(t) \lambda
$$
</div>

Where clade N' is the clade made up of the combination of two sister clades N and M.

For the parent clade to have gone extinct before the present, both daughters must have gone extinct, so:

(eq. 11.17)
<div>
$$
E_{N'}(t) = E_{N}(t) E_{M}(t)
$$
</div>

To apply this approach across an entire phylogenetic tree, we multiply equations 11.15 and 11.16 across all branches and nodes in the tree. Thus, the full likelihood is:

(eq. 11.18)
<div>
$$
L = \lambda^n[\prod_{k = 1}^{2N} e^{(\lambda-\mu)(t_{k,b} - t_{k,t})} \cdot \frac{(\lambda - (\lambda-\mu)e^{(\lambda - \mu)t_N})^2}{(\lambda - (\lambda-\mu)e^{(\lambda - \mu)t})^2}]
$$
</div>

Most methods fitting birth-death models to trees condition on the existence of a tree, which requires that the two lineages following the initial split in the tree survived to the present day. Conditioning on this, the likelihood becomes:


(eq. 11.19)
<div>
$$
L = (N-1)! \lambda^(N-2) \Big(\prod_{i=2}^{N-1} P_s(t_i, T)) \cdot (1-u_s(t_1))^2 \prod_{i=2}^{N-1} (1-u_s(t_i))
$$
</div>

Here $N$ is the number of tips in the tree and $t_i$ the speciation times as defined above. $P_s(t_i, T)$ and $u_s(t_i)$ are functions. $P_s(t_i, T)$ is the probability that a lineage at time $t_i$ leaves at least one lineage at the present:

(eq. 11.20)
<div>
$$
P_s(t_i, T) =
$$
</div>


and $u_s(t_i)$ is:

(eq. 11.21)
<div>
$$
u_s(t_i) =
$$
</div>



### Section 11.4b: Using maximum likelihood to fit a birth-death model

Given equation 11.11 for the likelihood, we can estimate birth and death rates using both ML and Bayesian approaches. For the ML estimate, we maximize equation 11.11 over $\lambda$ and $\mu$. For a pure-birth model, we can set $\mu$ = 0, and the maximum likelihood estimate of $\lambda$ can be calculated analytically as:

(eq. 11.14)

<div>
$$
\lambda=  \frac{n-2}{s_{branch}}
$$
</div>

where $s_{branch}$ is the sum of branch lengths in the tree,

(eq. 11.15)

<div>
$$
s_{branch} = \sum_{i=1}^{n-1} i t_i
$$
</div>


Equation 11.14 is also called the Kendal-Moran estimator of the speciation rate [@Nee2006-cd].

For a birth-death model, we can use numerical methods to maximize the likelihood over $\lambda$ and $\mu$.

For example, we can use ML to fit a birth-death model to the *Lupinus* tree [@Drummond2012-zs], which has 137 tip species and a total age of 16.6 million years. Doing so, we obtain ML parameter estimates of $\lambda = 0.46$ and $\mu = 0.20$, with a log-likelihood of $lnL_{bd} = 262.3$. Compare this to a pure birth model on the same tree, which gives $\lambda = 0.35$ and $lnL_{pb} = 260.4$. One can compare the fit of these two models using AIC scores: $AIC_{bd} = -520.6$ and $AIC_{pb} = -518.8$, so the birth-death model has a better (lower) AIC score but by less than two AIC units. A likelihood ratio test, which gives $\Delta = 3.7$ and $P = 0.054$. In other words, we estimate a non-zero extinction rate in the clade, but the evidence supporting that model over a pure-birth is not particularly strong.


### Section 11.4c: Using Bayesian MCMC to fit a birth-death model

We can also estimate birth and death rates using a Bayesian MCMC. We can use exactly the method spelled out above for clade ages and diversities, but substitute equation 11.11 for the likelihood, thus using the waiting times derived from a phylogenetic tree to estimate model parameters.

Applying this to Lupines with the same priors as before, we obtain the posterior distributions shown in figure 11.4. The mean of the posterior for each parameter is $\lambda = 0.48$ and $\mu = 0.23$, quite close to the ML estimates for these parameters.

![Figure 11.4. Posterior distribution for $b$ and $d$ for *Lupinus*  [@Drummond2012-zs].]({{ site.baseurl }}/images/figure11-4.png)

## Section 11.5: Sampling and birth-death models

It is important to think about sampling when fitting birth-death models to phylogenetic trees. If any species are missing from your phylogenetic tree, they will lead to biased parameter estimates. This is because missing species are disproportionally likely to connect to the tree on short, rather than long, branches. If we randomly sample lineages from a tree, we will end up badly underestimating both speciation and extinction rates (and wrongly inferring slowdowns; see chapter 12).

Fortunately, the mathematics for incomplete sampling of reconstructed phylogenetic trees has also been worked out. We can substitute in equations that include $\rho$, the proportion of sampled species, for equations 11.12 and 11.13 above:

(eq. 11.15)

<div>
$$
p_0 (t) = 1 - \frac{\rho (\lambda - \mu)}{\rho \lambda + (\lambda (1 - \rho)-\mu)e^{-(\lambda - \mu)t}}
$$
</div>

(eq. 11.16)

<div>
$$
p_1 (t) = \frac{\rho (\lambda - \mu)^2 e^{-(\lambda - \mu)t}}{(\rho \lambda + (\lambda (1 - \rho)-\mu)e^{-(\lambda - \mu)t})^2}
$$
</div>

##  Section 11.6: Summary

In this chapter, I described how to estimate parameters from birth-death models using data on species diversity and ages, and how to use patterns of tree balance to test hypotheses about changing birth and death rates. I also described how to calculate the likelihood for birth-death models on trees, which leads directly to both ML and Bayesian methods for estimating birth and death rates. Next, we will explore elaborations on birth-death models, and discuss models that go beyond constant-rates birth-death models to analyze the diversity of life on Earth.

## References
