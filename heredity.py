import csv
import itertools
import sys

PROBS = {

    # Unconditional probabilities for having gene
    "gene": {
        2: 0.01,
        1: 0.03,
        0: 0.96
    },

    "trait": {

        # Probability of trait given two copies of gene
        2: {
            True: 0.65,
            False: 0.35
        },

        # Probability of trait given one copy of gene
        1: {
            True: 0.56,
            False: 0.44
        },

        # Probability of trait given no gene
        0: {
            True: 0.01,
            False: 0.99
        }
    },

    # Mutation probability
    "mutation": 0.01
}


def main():

    # Check for proper usage
    if len(sys.argv) != 2:
        sys.exit("Usage: python heredity.py data.csv")
    people = load_data(sys.argv[1])

    # Keep track of gene and trait probabilities for each person
    probabilities = {
        person: {
            "gene": {
                2: 0,
                1: 0,
                0: 0
            },
            "trait": {
                True: 0,
                False: 0
            }
        }
        for person in people
    }

    # Loop over all sets of people who might have the trait
    names = set(people)
    for have_trait in powerset(names):

        # Check if current set of people violates known information
        fails_evidence = any(
            (people[person]["trait"] is not None and
             people[person]["trait"] != (person in have_trait))
            for person in names
        )
        if fails_evidence:
            continue

        # Loop over all sets of people who might have the gene
        for one_gene in powerset(names):
            for two_genes in powerset(names - one_gene):

                # Update probabilities with new joint probability
                p = joint_probability(people, one_gene, two_genes, have_trait)
                update(probabilities, one_gene, two_genes, have_trait, p)

    # Ensure probabilities sum to 1
    normalize(probabilities)

    # Print results
    for person in people:
        print(f"{person}:")
        for field in probabilities[person]:
            print(f"  {field.capitalize()}:")
            for value in probabilities[person][field]:
                p = probabilities[person][field][value]
                print(f"    {value}: {p:.4f}")


def load_data(filename):
    """
    Load gene and trait data from a file into a dictionary.
    File assumed to be a CSV containing fields name, mother, father, trait.
    mother, father must both be blank, or both be valid names in the CSV.
    trait should be 0 or 1 if trait is known, blank otherwise.
    """
    data = dict()
    with open(filename) as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row["name"]
            data[name] = {
                "name": name,
                "mother": row["mother"] or None,
                "father": row["father"] or None,
                "trait": (True if row["trait"] == "1" else
                          False if row["trait"] == "0" else None)
            }
    return data


def powerset(s):
    """
    Return a list of all possible subsets of set s.
    """
    s = list(s)
    return [
        set(s) for s in itertools.chain.from_iterable(
            itertools.combinations(s, r) for r in range(len(s) + 1)
        )
    ]


def joint_probability(people, one_gene, two_genes, have_trait):
    """
    Compute and return a joint probability.

    The probability returned should be the probability that
        * everyone in set `one_gene` has one copy of the gene, and
        * everyone in set `two_genes` has two copies of the gene, and
        * everyone not in `one_gene` or `two_gene` does not have the gene, and
        * everyone in set `have_trait` has the trait, and
        * everyone not in set` have_trait` does not have the trait.
    """
    #  P1 =  iterate over one_gene
    #       if person has no parent then no dependency
    #       use unconditional  probabilty given
    #        if person has parent then compute probabliy given parent has 1 gene, 2 gene , or 0 gene
    #       if has parents
    #          if parent has 2 genes then 1 gene
    #          if parent has 1 gene then prob 1/2
    #          if parent has 0 genes then 0 genes
    #       add for both parents consider joint prob since  combinqation requireed
    #       also for each scneario include mutation

    #      for have traint multiply the gene prob with the conditiounal trait prob


    # for now assume person has no parent or both parents.
    p_1_gene = 1
    for name in one_gene:
    p_1_gene = p_1_gene * get_gene_prob(people, name, 1)

    p_2_gene = 1
    for name in two_genes:
    p_2_gene = p_2_gene * get_gene_prob(people, name, 2)

    p_0_gene = 1
    for name in people:
        if name not in one_gene and name not in two_gene
            p_0_gene = p_0_gene * get_gene_prob(people, name, 0)

    p_trait = 1

    for name in have_trait
        person = people[name]
        p_trait = p_trait * (
            get_gene_prob(people, name, 0) * PROBS["trait"][0]["True"]
            + get_gene_prob(people, name, 1) * PROBS["trait"][1]["True"]
            + get_gene_prob(people, name, 2) * PROBS["trait"][2]["True"]


    p_no_trait = 1
    for name in have_trait
        person = people[name]
        p_trait = p_trait * (
            get_gene_prob(people, name, 0) * PROBS["trait"][0]["False"]
            + get_gene_prob(people, name, 1) * PROBS["trait"][1]["False"]
            + get_gene_prob(people, name, 2) * PROBS["trait"][2]["False"]

    return p_1_gene * p_2_gene * p_0_gene * p_trait * p_no_trait


def get_gene_prob(people, name, no_genes)
    Person person =  people[name]
    mother = person["mother"]
    father = person["father"]

    # this has to be recursive . we may have to go up generations to reach
    # the unconditional prob
    if mother is None and father is None:
        p_1_gene = PROBS["gene"][no_genes]    # unconditional prob
    else:
        # write all combinations
        #   Child Count M_Count     F_Count     M_Contrib    F_Contrib
        switch no_genes:
        case 0:
            ( get_gene_prob(people, mother, 0) * (1 - PROBS["mutation"]) +
              get_gene_prob(people, mother, 1) * 0.5 +
              get_gene_prob(people, mother, 2) * PROBS["mutation"]) )  *
            ( get_gene_prob(people, father, 0) * (1 - PROBS["mutation"]) +
              get_gene_prob(people, father, 1) * 0.5 +
              get_gene_prob(people, father, 2) * PROBS["mutation"]  )
        case 1:
            # try to make this easy
            # Prob mother is 1 & Father is 0
            # +  mother is 0 & father is 1
            # in each case start with 0 genes and go up to 2. consider mutation also. so 2 items for each
            ( get_gene_prob(people, mother, 0) * PROBS["mutation"] +
              get_gene_prob(people, mother, 1) * 0.5 +
              get_gene_prob(people, mother, 2) * (1- PROBS["mutation"])  )  *
            ( get_gene_prob(people, father, 0) * (1 - PROBS["mutation"]) +
              get_gene_prob(people, father, 1) * 0.5 +
              get_gene_prob(people, father, 2) * PROBS["mutation"]  )
        case 2:
             # +  mother is 1 & father is 1
            ( get_gene_prob(people, mother, 0) * PROBS["mutation"] +
              get_gene_prob(people, mother, 1) * 0.5 +
              get_gene_prob(people, mother, 2) * (1- PROBS["mutation"])  )  *
            ( get_gene_prob(people, father, 0) * PROBS["mutation"] +
              get_gene_prob(people, father, 1) * 0.5 +
              get_gene_prob(people, father, 2) *  (1- PROBS["mutation"])  )


            # if mother is  None and father is None:
            #     p1_gene = get_gene_prob(people, mother, 0) * PROBS["mutation"] * get_gene_prob(people, father, 0) * (1 - PROBS["mutation"]) +
            #                  get_gene_prob(people, mother, 0) * (1 - PROBS["mutation"]) * get_gene_prob(people, father, 0) * PROBS["mutation"]
            # if mother is not None and father is None:
            #     p1_gene = get_gene_prob(people, mother, 2) * (1- PROBS["mutation"]) *
            #                get_gene_prob(people, father, 0) * PROBS["mutation"]
            #                + get_gene_prob(people, mother, 1) * 0.5  # father cancels out
            # if mother is None and father is not None:
            #     p1_gene = (get_gene_prob(people, father, 2) * (1- PROBS["mutation"]) *
            #                 get_gene_prob(people, mother, 0) * PROBS["mutation"]
            #                + get_gene_prob(people, father, 1) * 0.5 # mother cancels out
            # if mother is not None and father is not None:
            #     p1_gene = (get_gene_prob(people, mother, 2) * (1- PROBS["mutation"])
            #                 * get_gene_prob(people, father, 1) * 0.5) # mutation prob cancels out for the father. not for the mother
            #                 + (get_gene_prob(people, mother, 2) * (1- PROBS["mutation"])

            #             + get_gene_prob(people, mother, 1) *   get_gene_prob(people, father, 1) # not the mutation prob cancels out
            #             +  (get_gene_prob(people, father, 2) * (1- PROBS["mutation"])
            #                 * get_gene_prob(people, mother, 1) * 0.5) # mutation prob cancels out for other way around here
            #               +  (get_gene_prob(people, father, 2) * (1- PROBS["mutation"])
            #                 * get_gene_prob(people, mother, 0)  *  (PROBS["mutation"])  # don't mutation for both father & mother
            # # need to calculate combination



def update(probabilities, one_gene, two_genes, have_trait, p):
    """
    Add to `probabilities` a new joint probability `p`.
    Each person should have their "gene" and "trait" distributions updated.
    Which value for each distribution is updated depends on whether
    the person is in `have_gene` and `have_trait`, respectively.
    """
    for person in one_gene:
        


def normalize(probabilities):
    """
    Update `probabilities` such that each probability distribution
    is normalized (i.e., sums to 1, with relative proportions the same).
    """
    raise NotImplementedError


if __name__ == "__main__":
    main()
