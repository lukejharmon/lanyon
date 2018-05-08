# Chapter 13: Characters and diversification rates

[pdf version]({{ site.baseurl }}/pdf/chapter12_beyondbd.pdf)

NOTE THIS CHAPTER IS IN PROGRESS!

## Section 12.1: The evolution of self-incompatibility

Most people have not spent a lot of time thinking about the sex lives of plants. The classic mode of sexual reproduction in angiosperms (flowering plants) involves pollen (the male gametophyte stage of the plant life cycle). Pollen lands on the pistil (the female reproductive structure) and produces a pollen tube. Sperm cells move down the pollen tube, and one sperm cell unites with the egg to form a new zygote in the ovule.

As you might imagine, plants have little control over what pollen grains land on their pistil (although plant species do have some remarkable adaptations to control pollination by animals; see xxx). In particular, this mode of reproduction leaves open the possibility of self-pollination, where pollen from a plant fertilizes eggs from the same plant. Self-fertilization (sometimes called selfing) is a form of asexual reproduction, but one that involves meiosis; as such, there are costs to self-fertilization. The main cost is inbreeding depression, a reduction in offspring fitness associated with recessive deleterious alleles across the genome.

Some species of angiosperms can avoid self-fertilization through self-incompatibility. In plants with self-incompatibility, the process by which the sperm meets the egg is interrupted at some stage if pollen grains have a genotype that is the same as the parent. This prevents selfing – and also prevents sexual reproduction with plants that have the same genotype(s) at loci involved in the process.

Species of angiosperms are about evenly divided between these two states Furthermore, self-incompatible species are scattered throughout the phylogenetic tree of angiosperms.

The evolution of selfing is a good example of a trait that might have a strong effect on diversification rates. One can easily imagine, for example, how incompatibility loci might facilitate the evolution of reproductive isolation among populations, and how lineages with such loci might diversity at a very different tempo than those without.

In this chapter, we will learn about a family of models where traits can affect diversification rates. I will also address some of the controversial aspects of these models and how we can improve these approaches in the future.

Key Biological Questions

	How does the behavior of birth-death models change when there are some character states that are associated with different rates of speciation and/or extinction compared to others?
	How we detect character-dependent diversification using comparative data?
	When and why might character-dependent approaches lead us astray?

## Section 12.2: A State-Dependent Model of Diversification

The models that we will consider in this chapter include both a model of trait evolution and an associated model of lineage diversification. In the simplest case, we can consider a model where the character has two states, 0 and 1. We need to model the transitions among these states, which we can do in an identical way to what we did in chapter xxx using a continuous-time Markov model. We express this model using two rate parameters, a forward rate q01 and a backwards rate q10.

We now consider the idea that diversification rates might depend on the character state. We assume that species with character state 0 have a certain speciation rate (lambda0) and extinction rate (mu0), and that species in 1 have potentially different rates of both speciation (lambda1) and extinction (mu1). Thus, we have a six-parameter model. We assume that parent lineages give birth to daughters with the same character state, that is that character states do not change at speciation.

It is straightforward to simulate evolution under our state-dependent model of diversification. We proceed in the same way as we did for birth-death models, by drawing waiting times, but these waiting times can be waiting times to the next character state change, speciation, or extinction event. In particular, imagine that there are n lineages present at time t, and that k of these lineages are in state 0 (and n-k are in state 1). The waiting time to the next event will follow an exponential distribution with a rate parameter of:

ρ=n∙(q_01+λ_0+μ_0 )+(k-n)∙(q_10+λ_1+μ_1 )

Once we have a waiting time, we can assign an event type depending on probabilities. For example, the probability that the event is a character state change from 0 to 1 is:

p_q01=(n∙q_01)/ρ

And the probability that the event is the extinction of a lineage with character state 1 is:

p_μ1=((n-k)∙μ_1)/ρ


Once we have picked an event in this way, we can randomly assign it to one of the lineages in the appropriate set, with each lineage equally likely to be chosen. We then proceed forwards in time until we have a dataset with the desired size or total time depth.

An example simulation is shown in Figure 14.xxx. As you can see, under these model parameters the impact of character states on diversification is readily apparent. In the next section we will figure out how to extract that information from our data.

Calculating Likelihoods for State-dependent diversification models

To calculate likelihoods for state-dependent diversification models we use a pruning algorithm with calculations that progress back through the tree from the tips to the root. The description of this algorithm in Maddison et al. (xxx) represents one of the clearest explanations of calculating likelihoods for comparative methods, and is a general approach that could be applied to a wide range of other models; thus, it is worth describing in some detail.

We consider a phylogenetic tree with data on the character states at the tips. For the purposes of this example, we will assume that the tree is complete and correct – we are not missing any species, and there is no phylogenetic uncertainty. We will come back to these two assumptions later in the chapter.

We need to obtain the probability of obtaining the data given the model (the likelihood). As we have seen before, we will calculate that likelihood going backwards in time using a pruning algorithm. The key principle here is that if we know the probabilities at some point in time on the tree, we can calculate those probabilities at some time point immediately before. By applying this method successively, we can move back towards the root of the tree. We move backwards along each branch in the tree, merging these calculations at nodes.
When we get to the root, we have the probability of the data given the model and the entire tree – that is, we have the likelihood. The other essential piece is that we have a starting point. When we start at the tips of the tree, we assume that our character states are fixed and known. We use the fact that we know all of the species and their character states at the present day as our starting point, and move backwards from there.For example, for a species with character state 0, the likelihood for state 0 is 1, and for state 1 is zero. In other words, at the tips of the tree we can start our calculations with a probability of 1 for the state that matches the tip state, and 0 otherwise.

This discussion also highlights the fact that incorporating uncertainty and/or variation in tip states for these algorithms is not computationally difficult – we just need to start from a different point at the tips. For example, if we are completely unsure about the tip state for a certain taxa, we can begin with likelihoods of 0.5 for starting in state 0 and 0.5 for starting in state 1. However, such calculations are not commonly implemented in comparative methods software.

We now need to consider the change in the likelihood as we step backwards through time in the tree. Instead of moving in discrete steps, as described in the toy example previously, we will consider some very small time interval ∆t, and use differential equations to find out what happens in the limit as this interval goes to zero (Figure 14.1). We will assume that the time interval ∆t is so small that, at most, one event (speciation or character change) has happened in that interval, but never more than one. We will consider the probability of the observed data given that the character is in each state at time t: p_0 (t) and p_1 (t). Assuming we know these probabilities, we calculate updated probabilities at some earlier time t+∆t: p_0 (t+∆t) and p_1 (t+∆t).



To calculate p_0 (t+∆t) and p_1 (t+∆t), we consider all of the possible things that could happen in a time interval ∆t along a branch in a phylogenetic tree that are compatible with our dataset (figure 14.2). First, nothing at all could have happened; second, our character state could have changed; and third, there could have been a speciation event. This last event might seem incorrect, as we are only considering changes along branches in the tree and not at nodes. If we did not reconstruct any speciation events at some point along a branch, then how could one have taken place? The answer is that a speciation event could have occurred but all taxa descended from that branch have since gone extinct. We must also consider the possibility that either the right or the left lineage went extinct following the speciation event; that is why the speciation event probabilities appear twice in Figure 14.2.



We can write an equation for these updated probabilities. We will consider the probability that the character is in state 0 at time t+∆t; the equation for state 1 is similar.

p_0 (t+∆t)=(1-μ_0 )∆t∙[(1-q_01 ∆t)(1-λ_0 ∆t) p_0 (t)+q_01 ∆t(1-λ_0 ∆t) p_1 (t)+2∙(1-q_01 ∆t) λ_0 ∆t∙E_0 (t)p_0 (t)]

We can multiply through and simplify. We will also drop any terms that include 〖∆t〗^2, which become vanishingly small as ∆t decreases. Doing that, we obtain:

p_0 (t+∆t)=[1-(λ_0+μ_0+q_01 )∆t] p_0 (t)+(q_01 ∆t) p_1 (t)+2(λ_0 ∆t) E_0 (t)p_0 (t)

p_1 (t+∆t)=[1-(λ_1+μ_1+q_10 )∆t] p_1 (t)+(q_10 ∆t) p_0 (t)+2(λ_1 ∆t) E_1 (t)p_1 (t)

We can then find the rate of change for these two equations by taking


We also need to know E_0 (t) and E_1 (t). These represent the probability that a lineage with state 0 or 1, respectively, and alive at time t will go extinct before the present day. Neglecting the derivation of these formulas, which can be found in Maddison et al. 2007, we have:

(dE_0)/dt=μ_0-(λ_0+μ_0+q_01 ) E_0 (t)+q_01 E_1 (t)+λ_0 〖E_0 (t)〗^2

(dE_1)/dt=μ_1-(λ_1+μ_1+q_10 ) E_1 (t)+q_10 E_0 (t)+λ_1 〖E_1 (t)〗^2




Along a single branch in a tree, we can sum together many such small time intervals. But what happens when we get to a node? Well, if we consider the time interval that countains the node, then we already know what happened – a speciation event. We also know that the two daughters immediately after the speciation event were identical in their traits (this is an assumption of the model). So we can calculate the likelihood for their ancestor for each state as the product of the likelihoods of the two daughter branches coming into that node. In this way, we merge our likelihood calculations along each branch when we get to tree nodes.

When we get to the root of the tree, we are almost done – but not quite! We have partial likelihood calculations for each character state – so we know, for example, the likelihood of the data if we had started with a root state of 0, and also if we had started at 1. To merge these we need to use probabilities of each character state at the root of the tree. For example, if we do not know the root state from any outside information, we might consider root probabilities for each state to be equal, 0.5 for state 0 and 0.5 for state 1. We then multiply the likelihood associated with each state with the root probability for that state. Finally, we add these likelihoods together to obtain the full likelihood of the data given the model.

The question of which root probabilities to use for this calculation has been discussed in the literature, and does matter in some applications. Other options include…

I have described the situation where we have two character states, but this method generalizes well to multi-state characters. We can describe the evolution of the character in the same way as described for multi-state discrete characters in chapter xxx. We then can assign unique diversification rate parameters to each of the k character states:. It is worth keeping in mind, though, that it is not too hard to construct a model where parameters are not identifiable and model fitting and estimation become very difficult.

## Section 12.3: ML and Bayesian Tests for State-Dependent Diversification

Now that we can calculate the likelihood for state-dependent diversification models, formulating ML and Bayesian tests follows the same pattern we have encounted before. For ML, some comparisons are nested and so you can use likelihood ratio tests. For example,

Potential Pitfalls and How to Avoid Them

Recently, a few papers have been published that are deeply critical of state-dependent diversification models. These papers raise substantive critiques that are critical to address when applying the methods described in this chapter to empirical data. In this section I will attempt to describe the critiques and their potential remedies.

The most serious limitation of state-dependent models as currently implemented is that they consider only a relatively small set of possible models. In particular, the approach we describe above compares two models: first, a model where birth and death rates are constant and do not depend on the state of the character; and second, a model where birth and death rates depend only on the character state. But there is another possibility that might be (in general) more common than either of the models we consider: birth and death rates vary, but in a way that is not dependent on the particular character we have chosen to analyze. I say that this is probably a common pattern because we know that birth and death rates vary tremendously across lineages in the tree of life, and it seems probable to me that many of our hypotheses about which characters might contribute to that variation are, at this point, stabs in the dark.

This issue is a normal one for statistical analyses – after all, there are always other models outside of our set of considered possibilities. However, in this case, the fact that state-dependent diversification models fail to consider the possibility outlined above causes a very particular – and peculiar – problem: if we apply the tests to empirical phylogenetic trees, even with made-up data, we almost always find statistically significant results. For example, Rabosky and Goldberg found that there is very often a statistically significant “signal” that the number of letters in a species name is significantly associated with speciation rates across a range of empirical datasets. This result might seem ridiculous and puzzling, as there is no way that species name length is associated with diversification processes. However, if we return to our alternative model above, then the results make sense. Rabosky and Goldberg simulated character evolution on real phylogenetic trees, and their results do not hold when the trees are simulated along with the characters (this is also why Rabosky and Goldberg’s results do not represent “type I errors,” contra their paper, because the data are not simulated under the null hypothesis). On these real trees, speciation and/or extinction rates have vary across clades. Among the two models that the authors consider, both are wrong; speciation and extinction are independent of the character but not constant through time. Of the two alternatives, the state-dependent model tends to fit better because, from a statistical point of view, it is important for the model to capture some variation in birth and death rates across clades. Even a random character will pick up some of this variation, so that the alternative model tends to fit better than the null – even though, in this case, the character has nothing to do with diversification!

Fortunately, there are a number of ways to deal with this problem. First, one can compare the statistical support for the state-dependent model with the support that one obtains for random data. The random data could be simulated on the tree, or one could permute the tips or draw random data from a multinomial distribution. One can then compare, for example, the distribution of dAICc scores obtained from these permutations to the dAICc for the original data. Alternatively, we could explicitly consider the possibility that some unmeasured character is actually the thing that is influencing diversification rates. This latter approach is the most elegant as we can directly add the model described in this section to our list of candidates (see HISSE paper for more details).

A more general critique of state-dependent models of diversification was raised by Maddison and Fitzjohn (Xxx). This paper pointed out that statistically significant results for these tests can be driven by an event on a single branch of a tree, and therefore be unreplicated. This is a good criticism that applies equally well to a range of comparative methods. I will address this criticism later in the book when I discuss model adequacy for comparative methods.

## Section 12.4: Summary

Many evolutionary models postulate a link between species characteristics and speciation, extinction, or both. These hypotheses can be tested using state-dependent diversification models, which explicitly consider the possibility that species’ characters affect their diversification rates. State-dependent models as currently implemented have some potential problems, but there are methods to deal with these critiques. The overall ability of state-dependent models to explain broad patterns of evolutionary change remains to be determined, but represents a promising avenue for future research.


## Footnotes


# References
