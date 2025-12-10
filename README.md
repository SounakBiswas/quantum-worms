Attempts to construct quantum worms.

### What is the myopic-equivalent here.

We have two kinds of vertices, typ-0 has 3 and 1 dimers with equal weight. typ-1 has only 1-dimer. And then we have another kind of hyperedge vertices. which has 0 or 1 dimer touching it, call it typ-2.  In the following D(v1,v2)=0/1 denotes that the edge (v1,v2) has or does not have a dimer.

1. Start at a random entry vertex n. If n is typ-1 then then must enter the next vertex v such that D(n,v)=1 . If n is typ-0, then D(n,v)=1 or 0, choose whatever, randomly.
    1. If D(n,v)=1: v is typ-1 then x can be chosen uniformly among D(v,x)=0 links. If v is typ-0, Then D(v,x)=1 necessarily.
    2. If D(n,v)=0: v is typ-1 and x is chosen uniquely in D(v,x)=1. If v is typ-0 then choose any x
2. If D(v,x)=0, x is typ-1 then v1 is chosen such that D(x,n1)=1. anything if x is typ-0
3. If D(v,x)=1 x is typ-1 then v1 is chosen randomly. 


There is no need for myopic. All steps respect det-bal



### What to do about the hyper edges?

Don't touch the hyper-edges


### What are the main issues still to be addressed

1. How to reconstruct spins from dimer configs
2. Stitch together dual graph accross the imaginary-time boundary
3. What is the best way to backtrack loops that violate winding number parity
