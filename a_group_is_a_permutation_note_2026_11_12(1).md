# A group is a permutation: the usual representation theory from mirrors, turns and blocks

*What this is.* A separate note, not part of the master and not a draft of a section. It takes one way of playing with groups and follows it on small matrices, written out. It ends at the standard facts of representation theory, each read off a matrix. Sections 1 to 8 follow the game on examples; sections 9 to 11 turn it into an algorithm, derive the usual formulas from it, and run it. Section 12 lets the angle be anything: the ring becomes a circle, the whole numbers or a line, and the rotations of space are done the same way. Every statement is checked in `scripts/groups_as_permutations.py` (9 checks) and `scripts/group_blocks_algorithm.py` (the algorithm, run on seven groups) and `scripts/fixed_angle_to_any_angle.py` (14 checks, section 12).

**The game, in four moves.**

1. A group is a permutation.
2. A permutation is made of swaps, and a swap is a mirror.
3. In one more dimension a mirror is a half-turn, so everything becomes a turn.
4. Turns are made diagonal by Fourier. Do that, then come back down.

Moves 1 to 3 are right as they stand. Move 4 is right when all the operations commute. When they do not, the most that can be done is *blocks*, and the sizes of the blocks are the content of the subject. That one correction is what turns the game into the usual theory.

## 1. A group is a permutation

*The picture first.* Three things in a ring. The only things to do are: leave them, shift them one place, shift them two places. Call these $e$, $r$, $r^2$. Shifting moves thing 1 to place 2, thing 2 to place 3, thing 3 to place 1. As a matrix on the three places:

$$r=\begin{pmatrix}0&0&1\\1&0&0\\0&1&0\end{pmatrix},\qquad r^2=r\,r,\qquad r^3=e .$$

1. one entry $1$ in each row and each column: a permutation.

Any group can be written this way, by letting it act on its own elements. Section 5 does it for a group of six.

## 2. A permutation is made of swaps, and a swap is a mirror

The swap of things 1 and 2, and the swap of things 2 and 3:

$$S_{12}=\begin{pmatrix}0&1&0\\1&0&0\\0&0&1\end{pmatrix},\qquad S_{23}=\begin{pmatrix}1&0&0\\0&0&1\\0&1&0\end{pmatrix},\qquad S_{12}\,S_{23}=\begin{pmatrix}0&0&1\\1&0&0\\0&1&0\end{pmatrix}=r .$$

2. the shift is two swaps. A swap done twice is nothing, $S^2=1$, and its determinant is $-1$: it is a mirror. $S_{12}$ is the mirror in the plane $x_1=x_2$: it keeps every point of that plane and reverses the direction $(1,-1,0)$ across it.

3. two mirrors, one after the other, are a turn. The shift $r$ keeps the diagonal $(1,1,1)$ and turns everything else about it by a third of a turn, $120^\circ$. (The two mirror planes $x_1=x_2$ and $x_2=x_3$ meet at $60^\circ$ along the diagonal; two mirrors at an angle make a turn by twice the angle.) So "shift three things" and "turn by a third about the diagonal" are the same matrix.

## 3. One more dimension turns a mirror into a half-turn

Two things, and their swap. It is a mirror, determinant $-1$. Add one line, and let the swap reverse it:

$$S=\begin{pmatrix}0&1\\1&0\end{pmatrix}\qquad\longrightarrow\qquad\begin{pmatrix}0&1&0\\1&0&0\\0&0&-1\end{pmatrix} .$$

4. determinant $+1$. It keeps the direction $(1,1,0)$ and reverses the two directions across it: a half-turn about $(1,1,0)$. Seen from the plane of the first two coordinates it is still the mirror.

5. the same for any permutation $g$: put $\det g$ on the added line, $g\to\mathrm{diag}(g,\det g)$. Products are kept, and every element now has determinant $+1$. What the added line holds is $+1$ for a permutation made of an even number of swaps and $-1$ for an odd number: the *sign*. It comes back in section 6.

What this move is for: a turn in a plane, $\begin{pmatrix}\cos\theta&-\sin\theta\\\sin\theta&\cos\theta\end{pmatrix}$, becomes diagonal with the two numbers $e^{+i\theta}$ and $e^{-i\theta}$. A mirror has no angle to do that with.

## 4. Operations that commute: Fourier makes them all diagonal

*The picture first.* $n$ things in a ring. The group is the $n$ shifts $r^m$. They all commute. A pattern on the ring that repeats after a whole number of waves is only multiplied by a number when it is shifted. There are $n$ such patterns, one for each number of waves $k=0,\dots,n-1$.

With $\zeta=e^{2\pi i/n}$, the pattern with $k$ waves and what the shift does to it:

$$f_k=\frac1{\sqrt n}\begin{pmatrix}1\\\zeta^{k}\\\zeta^{2k}\\\vdots\\\zeta^{(n-1)k}\end{pmatrix},\qquad r\,f_k=\zeta^{-k}\,f_k,\qquad r^m\,f_k=\zeta^{-km}\,f_k .$$

6. in the basis of the $f_k$ every shift is diagonal. For three things, with $\zeta=e^{2\pi i/3}$:

$$r\ \longrightarrow\ \begin{pmatrix}1&0&0\\0&\zeta^{-1}&0\\0&0&\zeta^{-2}\end{pmatrix},\qquad r^2\ \longrightarrow\ \begin{pmatrix}1&0&0\\0&\zeta^{-2}&0\\0&0&\zeta^{-1}\end{pmatrix}\quad(\zeta^3=1).$$

7. $n$ blocks of size one. Block $k$ assigns to each element $r^m$ the number $\zeta^{-km}$. Writing any pattern as a sum of the $f_k$ is the Fourier series. For a group whose operations all commute, this is the whole of the theory: the blocks are the frequencies.

## 5. Operations that do not commute: blocks

*The picture first.* A triangle with equal sides. It can be turned by a third, by two thirds, or mirrored in one of three lines. With "leave it" that is six operations. Turn then mirror is not the same as mirror then turn. Two diagonal matrices always commute, so no basis can make all six diagonal. What can be done is to find a basis in which each of the six is made of small square blocks along the diagonal, the same sizes for all six.

*The six.* $r$ is the turn by a third, $s$ one of the mirrors. The six are $e$, $r$, $r^2$, $s$, $sr$, $sr^2$. Their table of products, row times column:

| left $\backslash$ right | $e$ | $r$ | $r^2$ | $s$ | $sr$ | $sr^2$ |
|---|---|---|---|---|---|---|
| $e$ | $e$ | $r$ | $r^2$ | $s$ | $sr$ | $sr^2$ |
| $r$ | $r$ | $r^2$ | $e$ | $sr^2$ | $s$ | $sr$ |
| $r^2$ | $r^2$ | $e$ | $r$ | $sr$ | $sr^2$ | $s$ |
| $s$ | $s$ | $sr$ | $sr^2$ | $e$ | $r$ | $r^2$ |
| $sr$ | $sr$ | $sr^2$ | $s$ | $r^2$ | $e$ | $r$ |
| $sr^2$ | $sr^2$ | $s$ | $sr$ | $r$ | $r^2$ | $e$ |


8. each row is a rearrangement of the six: the group is a permutation of its own elements. And $r\,s=s\,r^2$, which is not $s\,r$.

*As matrices.* Multiplying all six on the left by $r$, or by $s$, read off the rows $r$ and $s$ of the table:

$$\text{left }r=\begin{pmatrix}0&0&1&0&0&0\\1&0&0&0&0&0\\0&1&0&0&0&0\\0&0&0&0&1&0\\0&0&0&0&0&1\\0&0&0&1&0&0\end{pmatrix},\qquad\text{left }s=\begin{pmatrix}0&0&0&1&0&0\\0&0&0&0&1&0\\0&0&0&0&0&1\\1&0&0&0&0&0\\0&1&0&0&0&0\\0&0&1&0&0&0\end{pmatrix}.$$

9. two permutations of six things that do not commute.

*What does commute with everything.* Two kinds of thing commute with all six left multiplications.

- Every multiplication on the right, because $(g\,x)\,h=g\,(x\,h)$.
- The sum of the three mirrors, and the sum of the two turns. Renaming the corners of the triangle turns one mirror into another, so the sum of all three is the same whoever names the corners.

These can be made diagonal. The sum of the three mirrors takes the value $3$ on one direction, $-3$ on one direction, and $0$ on a space of four. That is the Fourier step for a group that does not commute: make diagonal the things that do commute. It sorts the six directions into three kinds, of sizes $1$, $1$ and $4$.

*The blocks.* The triangle itself gives six $2\times2$ matrices, the ones that move its corners (corners at $90^\circ$, $210^\circ$, $330^\circ$):

$$r:\ \begin{pmatrix}-\tfrac12&-\tfrac{\sqrt3}2\\\tfrac{\sqrt3}2&-\tfrac12\end{pmatrix},\qquad s:\ \begin{pmatrix}-1&0\\0&1\end{pmatrix},\qquad\text{and their products.}$$

Take as new basis of the six directions:

- direction 1: the same amount on all six elements;
- direction 2: $+1$ on the three turns $e,r,r^2$ and $-1$ on the three mirrors;
- directions 3 to 6: the four entries of the triangle's $2\times2$ matrix, each entry listed over the six elements.

In that basis:

$$\text{left }r\ \longrightarrow\ \begin{pmatrix}1&0&0&0&0&0\\0&1&0&0&0&0\\0&0&-\tfrac12&-\tfrac{\sqrt3}2&0&0\\0&0&\tfrac{\sqrt3}2&-\tfrac12&0&0\\0&0&0&0&-\tfrac12&-\tfrac{\sqrt3}2\\0&0&0&0&\tfrac{\sqrt3}2&-\tfrac12\end{pmatrix},$$

$$\text{left }s\ \longrightarrow\ \begin{pmatrix}1&0&0&0&0&0\\0&-1&0&0&0&0\\0&0&-1&0&0&0\\0&0&0&1&0&0\\0&0&0&0&-1&0\\0&0&0&0&0&1\end{pmatrix}.$$

10. all six become blocks of sizes $1$, $1$, $2$, $2$. The first block is always $1$. The second is $+1$ for a turn and $-1$ for a mirror: the sign of section 3. The third is the triangle's own $2\times2$ matrix, and the fourth is the same matrix again.

## 6. The usual facts, read off

11. **A block of size $d$ appears $d$ times.** The blocks of size one appear once each; the block of size two appears twice.

12. **The squares of the sizes add up to the size of the group**: $1^2+1^2+2^2=6$. The reason is on the right side. What commutes with all the left multiplications is exactly the right multiplications, and there are six of them. A thing that commutes with all the blocks can only mix identical blocks with each other, and mixing $d$ identical blocks takes a $d\times d$ array of numbers. So the count of right multiplications is $\sum d^2$.

13. **There are as many kinds of block as there are classes.** A class is a set of elements that renaming turns into each other. The triangle has three: $\{e\}$, the two turns, the three mirrors. And three kinds of block.

14. **The table of traces.** Take the trace of each block for each class:

| block | $e$ | the two turns | the three mirrors |
|---|---|---|---|
| one | $1$ | $1$ | $1$ |
| sign | $1$ | $1$ | $-1$ |
| triangle | $2$ | $-1$ | $0$ |

Any two rows are perpendicular when each class is counted with its size: for the last two, $1\cdot(1)(2)+2\cdot(1)(-1)+3\cdot(-1)(0)=0$. These rows do for the triangle what the frequencies do for the ring.

15. **The ring is the case where every block has size one**: $n=1+1+\dots+1$. Operations that commute have nothing to mix.

## 7. The same game on the walk's three mirrors

The three runs of the walk, $A_x$, $A_y$, $A_z$, are mirrors: each squares to $+1$, and any two anticommute. Their products are sixteen things, $\pm1$, $\pm A_k$, $\pm R$, $\pm\omega$.

16. two perpendicular mirrors are a half-turn: $R_{yz}=A_yA_z$, applied on both sides as $R_{yz}\,A\,R_{yz}^{-1}$, reverses $A_y$ and $A_z$ and keeps $A_x$. Done twice it is a full turn, and the product is $R_{yz}^2=-1$, not $+1$. The group made of the mirrors is twice as large as the group of turns they make, and $-1$ is the full turn.

17. the sixteen have ten classes, so ten kinds of block: eight of size one and two of size two, $16=8\times1^2+2\times2^2$. The blocks of size one do not tell $-1$ from $+1$. The part on which $-1$ acts as $-1$ has dimension eight: the two blocks of size two, each twice. That part is the eight detectors of the walk, and "each twice" is the original four and the added four. The details are in `detectors_as_a_group_draft_2026_11_12`.

## 8. What the game had wrong, and what it had right

- **"Make it diagonal" has to become "sort it into blocks".** Only operations that commute can all be diagonal together. What can always be made diagonal is what commutes with everything: the sums over classes. They tell the kinds of block apart, and the sizes of the blocks are what is left.
- **Raising the dimension does not change the group.** It makes every element a turn, which is what gives complex numbers something to make diagonal. The group is the same finite group.
- **"Coming back down" is taking off the sign.** The line added in section 3 carries exactly the block called *sign*. Nothing else was added.
- **What it had right is the route.** A group acting on itself, sorted into blocks, contains every block there is, each as many times as its size. That is the standard way to find them all.

## 9. The game as an algorithm

*What goes in.* Only the table of products of the group: $N$ elements, and for each pair which element their product is. *What comes out.* The sizes of the blocks, the blocks themselves as matrices, the table of traces, and the change of basis that sorts the group into blocks. Sections 5 and 6 did this by hand for the triangle; here is the same for any table.

*Symbols.* $e$ is the element that changes nothing. $e_y$ is the direction of element $y$ in a space of $N$ directions. A bar is the complex conjugate and $\dagger$ the conjugate transpose.

**Step 1. The group as permutations.** For each element $g$, two $N\times N$ matrices:

$$L(g)\,e_y=e_{g\,y},\qquad R(g)\,e_y=e_{y\,g^{-1}} .$$

18. $L(g)$ has a $1$ in row $x$, column $y$ when $x=g\,y$, and $0$ elsewhere; $R(g)$ has a $1$ when $x=y\,g^{-1}$. Both keep products, $L(g)L(h)=L(gh)$ and $R(g)R(h)=R(gh)$ (the $g^{-1}$ in $R$ is there for that), and every left commutes with every right.

**Step 2. The classes and their sums.** The class of $g$ is the set of all $h\,g\,h^{-1}$. Call the classes $K_1,\dots,K_m$, and add up the lefts of each:

$$C_j=\sum_{g\in K_j}L(g) .$$

19. each $C_j$ commutes with every left and every right, and the $C_j$ commute with each other.

**Step 3. Make the class sums diagonal, all together.** They commute with each other, so one basis makes them all diagonal. This is the only place where anything is "made diagonal". The $N$ directions fall into pieces $V_1,\dots,V_m$: on each piece every $C_j$ is a single number $\lambda_{ij}$, and two different pieces differ in at least one of these numbers. (In practice: take one combination $\sum_jt_j\,C_j$ with arbitrary numbers $t_j$ and find the spaces on which it is a number. The numbers are chosen so that the combination equals its own conjugate transpose, which makes it diagonal in a perpendicular basis: the class of the inverses gets the conjugate number.)

**Step 4. The size.** The piece $V_i$ has $d_i^{\,2}$ directions. Its square root $d_i$ is the size of the block of kind $i$.

**Step 5. Split a piece into copies.** Inside $V_i$ the rights still act, and they commute with the lefts. Take one combination of the rights with arbitrary numbers, $\sum_ht_h\,R(h)$ with $t_{h^{-1}}=\overline{t_h}$ (so that it equals its own conjugate transpose), and make it diagonal inside $V_i$. Its values come in $d_i$ sets of $d_i$ equal ones. The $d_i$ directions of one set are one copy of the block: every left multiplication sends them into combinations of themselves.

**Step 6. The block.** Write those $d_i$ directions, made perpendicular and of length one, as the columns of an $N\times d_i$ matrix $W_i$. The block of kind $i$ is

$$\rho_i(g)=W_i^{\dagger}\,L(g)\,W_i\qquad(d_i\times d_i),\qquad\rho_i(g)\,\rho_i(h)=\rho_i(gh),\qquad\chi_i(g)=\operatorname{Tr}\rho_i(g) .$$

**Step 7. The new basis.** For each kind $i$ and each entry $(a,b)$ of its block, one direction:

$$q_{i,ab}=\sqrt{\frac{d_i}{N}}\ \sum_g\overline{\rho_i(g)_{ab}}\ e_g .$$

20. $N$ directions in all, perpendicular and of length one. With them as the columns of $Q$,

$$Q^{\dagger}\,L(g)\,Q=\begin{pmatrix}\rho_1(g)&&&\\&\ddots&&\\&&\rho_i(g)&\\&&&\ddots\end{pmatrix},\qquad\text{each }\rho_i(g)\text{ repeated }d_i\text{ times.}$$

21. the arbitrary numbers of steps 3 and 5 only choose a basis inside a block; sizes and traces do not depend on them. The blocks come out complex in general. Over the real numbers a complex block and its conjugate join into one real block of twice the size: that is the walk's block of four counters, with $\omega$ as its quarter turn.

## 10. The usual formulas, each from the game

Two facts about a group acting on itself carry everything.

**(a) What commutes with every left is a combination of the rights.** Let $X$ commute with every $L(g)$, and write $X\,e_e=\sum_hx_h\,e_h$ for what it does to the direction of $e$. Then

$$X\,e_g=X\,L(g)\,e_e=L(g)\,X\,e_e=\sum_hx_h\,e_{g\,h}=\sum_hx_h\,R(h^{-1})\,e_g$$

$$\Longrightarrow\qquad X=\sum_hx_h\,R(h^{-1}) .$$

22. $N$ numbers $x_h$, and the $N$ rights are independent. In the same way, what commutes with every right is a combination of the lefts.

**(b) A left multiplication by anything but $e$ leaves no element in place.** $g\,y=y$ only for $g=e$. So

$$\operatorname{Tr}L(g)=\begin{cases}N,&g=e\\0,&g\neq e .\end{cases}$$

Write $n_i$ for the number of times the block of kind $i$ appears in the sorted form.

**10.1 The squares of the sizes add up to $N$, and a block of size $d$ appears $d$ times.** Three counts of the sorted form:

$$\sum_in_i\,d_i=N,\qquad\sum_in_i^{\,2}=N,\qquad\sum_id_i^{\,2}=N .$$

23. the first is the number of directions. The second: what commutes with the lefts can only mix identical blocks, an $n_i\times n_i$ array of numbers for each kind, and by (a) there are $N$ such things. The third: the $N$ lefts are independent, and between them they fill every block completely, a $d_i\times d_i$ array for each kind. (That a block which cannot be split further is filled completely is the one step taken here from the standard theory without proof.) Then

$$\sum_i(n_i-d_i)^2=\sum_in_i^{\,2}+\sum_id_i^{\,2}-2\sum_in_i\,d_i=N+N-2N=0,$$

and a sum of squares is zero only if each is: $n_i=d_i$.

**10.2 As many kinds of block as classes.** Count in two ways the things that commute with every left *and* every right. In the sorted form such a thing is one number for each kind of block. By (a) it is $\sum_hx_h\,R(h^{-1})$, and commuting with the rights as well makes $x_h$ the same on all of a class: one number for each class. So $m$ kinds and $m$ classes.

**10.3 The sizes weigh the traces to zero.** Take the trace of the sorted form and use (b) and $n_i=d_i$:

$$\sum_id_i\,\chi_i(g)=\begin{cases}N,&g=e\\0,&g\neq e .\end{cases}$$

**10.4 The traces are read off step 3.** $C_j$ commutes with everything, so on kind $i$ it is the number $\lambda_{ij}$ of step 3. Take the trace over one block of kind $i$:

$$\sum_{g\in K_j}\chi_i(g)=|K_j|\,\chi_i(K_j)=\lambda_{ij}\,d_i\qquad\Longrightarrow\qquad\chi_i(K_j)=\frac{d_i\,\lambda_{ij}}{|K_j|} .$$

**10.5 The projection on one kind.** Let $P_i$ be $1$ on the piece $V_i$ and $0$ on the others. It commutes with everything, so by (a) it is a combination of the lefts, $P_i=\sum_gc_i(g)\,L(g)$. Multiply by $L(h^{-1})$ and take the trace. On the left, $d_i$ copies of block $i$ give $d_i\,\chi_i(h^{-1})$. On the right, (b) leaves only $g=h$, which gives $N\,c_i(h)$. So

$$P_i=\frac{d_i}{N}\sum_g\chi_i(g^{-1})\,L(g) .$$

**10.6 The rows of the table of traces are perpendicular.** $P_i\,P_k=0$ for two different kinds and $P_i\,P_i=P_i$, whose trace is $d_i^{\,2}$. With 10.5 and (b), $\operatorname{Tr}(P_iP_k)=\frac{d_id_k}{N}\sum_g\chi_i(g^{-1})\,\chi_k(g)$. So

$$\frac1N\sum_g\chi_i(g^{-1})\,\chi_k(g)=\begin{cases}1,&i=k\\0,&i\neq k,\end{cases}\qquad\chi_i(g^{-1})=\overline{\chi_i(g)} .$$

24. a square table with perpendicular rows has perpendicular columns: $\sum_i\chi_i(g)\,\overline{\chi_i(h)}$ is $N/|K|$ when $g$ and $h$ are in the same class $K$, and $0$ otherwise.

**10.7 Fourier on a group.** A pattern is a number $f(g)$ on each element. Its coefficient on kind $i$ is a $d_i\times d_i$ matrix, and the pattern is recovered from the coefficients:

$$\hat f_i=\sum_gf(g)\,\rho_i(g),\qquad f(g)=\frac1N\sum_id_i\operatorname{Tr}\!\big(\rho_i(g^{-1})\,\hat f_i\big),$$

$$\sum_g|f(g)|^2=\frac1N\sum_id_i\operatorname{Tr}\!\big(\hat f_i^{\,\dagger}\hat f_i\big) .$$

25. these are step 7 read as a change of basis: $Q$ keeps lengths, and the component of the pattern on $q_{i,ab}$ is $\sqrt{d_i/N}$ times the entry $(a,b)$ of $\hat f_i$. And a combination of lefts becomes a product of coefficients:

$$\Big(\sum_xf(x)\,L(x)\Big)h\quad\text{has coefficients}\quad\hat f_i\,\hat h_i .$$

26. for $n$ things in a ring every $d_i=1$, $\rho_k(r^m)=\zeta^{-km}$, and these are the Fourier series of section 4 and its rule for products.

## 11. The algorithm run on seven groups

`scripts/group_blocks_algorithm.py` is steps 1 to 7 as a program. It is given only the table of products. For each group it finds the blocks and then checks 10.1 to 10.7 on what it found.

| group | elements | classes | blocks | check: $\sum d^2$ |
|---|---|---|---|---|
| six things in a ring | 6 | 6 | 6 of size 1 | $6\times1^2=6$ |
| the triangle (6 symmetries) | 6 | 3 | 2 of size 1, 1 of size 2 | $2\times1^2+2^2=6$ |
| the square (8 symmetries) | 8 | 5 | 4 of size 1, 1 of size 2 | $4\times1^2+2^2=8$ |
| the quaternion units (8) | 8 | 5 | 4 of size 1, 1 of size 2 | $4\times1^2+2^2=8$ |
| even permutations of four things (12) | 12 | 4 | 3 of size 1, 1 of size 3 | $3\times1^2+3^2=12$ |
| the walk's three mirrors (16) | 16 | 10 | 8 of size 1, 2 of size 2 | $8\times1^2+2\times2^2=16$ |
| all permutations of four things (24) | 24 | 5 | 2 of size 1, 1 of size 2, 2 of size 3 | $2\times1^2+2^2+2\times3^2=24$ |


27. the ring gives the six sixth roots of one, and the triangle gives the table of section 6. The square and the quaternion units are different groups with the same sizes of blocks and the same table of traces: the table of traces does not tell every group apart. The three mirrors of the walk give eight blocks of size one and two of size two, as in section 7.

## 12. From a fixed angle to any angle

*The picture first.* A ring of $N$ things is the turns by a fixed angle, a whole number of times $2\pi/N$. Let $N$ grow and the angle can be anything. What is held fixed while $N$ grows decides what the ring becomes: a circle, the whole numbers, or a line. In each case a sum over the group becomes an integral. What is lost in replacing a sum by an integral is not a small correction: it is the same integral again, one period further on, and again two periods on. These are *images*. With the images kept the passage is exact. The rotations of space then go through the same steps, and the finite count $\sum d^2=N$ becomes a count up to a cut-off.

The first half of this section is an older derivation of the Fourier integral from discrete steps, with its one wrong step replaced by the exact one.

*Symbols.* $N=2M+1$ is odd. The places are $n=-M,\dots,M$, at positions $t=n\epsilon$. The labels are $k=-M,\dots,M$, at frequencies $\omega=k\Delta$. The two steps are tied by $\epsilon\,\Delta\,N=2\pi$. $P=N\epsilon$ is the length of the ring, and $W'=N\Delta/2$ is half the width of the band of frequencies.

**(A) The ring, exactly.** The pattern with $k$ waves has the value $e^{ik\Delta\,n\epsilon}$ at place $n$ (section 4 in these names). Since $\epsilon\Delta/2\pi=1/N$,

$$\frac{\epsilon\Delta}{2\pi}\sum_{k=-M}^{M}e^{ik\Delta(m-n)\epsilon}=\begin{cases}1,&m=n\\0,&m\neq n\end{cases}\tag{1}$$

28. for $m\neq n$ the terms are $N$ points spaced evenly round a circle, and they add to zero. This is 10.3 for the ring, where every block has size one. So for any $N$ numbers $f(n\epsilon)$

$$f(m\epsilon)=\frac{\Delta}{2\pi}\sum_{k=-M}^{M}\Big[\epsilon\sum_{n=-M}^{M}f(n\epsilon)\,e^{-ik\Delta\,n\epsilon}\Big]e^{ik\Delta\,m\epsilon} .\tag{2}$$

29. exact, with no condition on the numbers.

**(B) The sum over the labels as a shift.** Write $\tau=(m-n)\epsilon$ and $D=d/d\omega$. The operator $e^{aD}$ shifts $\omega$ by $a$, so $e^{k\Delta D}e^{i\tau\omega}\big|_{\omega=0}=e^{ik\Delta\tau}$, and the sum over $k$ is a geometric series of shifts:

$$\sum_{k=-M}^{M}e^{k\Delta D}=e^{-M\Delta D}\,\frac{1-e^{N\Delta D}}{1-e^{\Delta D}}=\frac{e^{-N\Delta D/2}-e^{N\Delta D/2}}{e^{-\Delta D/2}-e^{\Delta D/2}}=\frac{\sinh(N\Delta D/2)}{\sinh(\Delta D/2)} .\tag{3}$$

30. the middle step multiplies the top and the bottom by $e^{-\Delta D/2}$ and uses $M+\tfrac12=N/2$. Acting on $e^{i\tau\omega}$ at $\omega=0$, $D$ becomes $i\tau$: the sum is $\sin(N\Delta\tau/2)/\sin(\Delta\tau/2)$.

**(C) An integral, and the rest.** Split (3) into two factors:

$$\Delta\sum_{k=-M}^{M}e^{k\Delta D}=\underbrace{\frac{2\sinh(W'D)}{D}}_{\text{the integral over }(-W',\,W')}\times\underbrace{\frac{\Delta D/2}{\sinh(\Delta D/2)}}_{\text{the rest}} .\tag{4}$$

31. the first factor on $e^{i\tau\omega}$ at $\omega=0$ is $2\sin(W'\tau)/\tau=\int_{-W'}^{W'}e^{i\omega\tau}\,d\omega$. The second is $(\Delta\tau/2)/\sin(\Delta\tau/2)$. It is close to one only while $\Delta\tau$ is small. But $\tau$ runs up to nearly the whole length $P=2\pi/\Delta$, where $\Delta\tau/2$ comes near $\pi$ and the factor grows without limit. It cannot be dropped, and it is not $1+O(\Delta)$. What it is, exactly, follows from

$$\frac1{\sin x}=\sum_{\ell=-\infty}^{\infty}\frac{(-1)^\ell}{x-\pi\ell},\qquad x=\frac{\Delta\tau}2,\qquad x-\pi\ell=\frac\Delta2\,(\tau-\ell P),$$

and from $\sin\!\big(W'(\tau-\ell P)\big)=\sin(W'\tau-\ell N\pi)=(-1)^\ell\sin(W'\tau)$ for odd $N$:

$$\boxed{\ \frac{\epsilon\Delta}{2\pi}\sum_{k=-M}^{M}e^{ik\Delta\tau}=\frac{\epsilon}{2\pi}\sum_{\ell=-\infty}^{\infty}\frac{2\sin\!\big(W'(\tau-\ell P)\big)}{\tau-\ell P}\ }\tag{5}$$

32. the term $\ell=0$ is the integral. The others are the same integral with $\tau$ moved by one length of the ring, two lengths, and so on: its images. A sum over labels is an integral over labels plus its images. With $N=25$, $\epsilon=0.37$ and $P=9.25$:

| $\tau$ | the exact sum | the integral's kernel alone | kernel plus images |
|---|---|---|---|
| $0.31$ | $+0.185627$ | $+0.185284$ | $+0.185627$ |
| $1.70$ | $+0.070071$ | $+0.066242$ | $+0.070071$ |
| $3.30$ | $+0.011194$ | $+0.008993$ | $+0.011194$ |
| $4.40$ | $-0.013365$ | $-0.008917$ | $-0.013365$ |
| $4.62$ | $+0.039964$ | $+0.025469$ | $+0.039964$ |


33. the same happens on the other side. If a pattern on a circle has coefficients $c_k$ for all whole numbers $k$, and it is sampled on a ring of $N$ places, the ring's coefficient $k$ is $\sum_\ell c_{k+\ell N}$: the coefficient plus its images, $N$ labels apart.

**(D) Three limits.** Put (5) into (2):

$$f(m\epsilon)=\epsilon\sum_{n=-M}^{M}f(n\epsilon)\ \frac1{2\pi}\sum_{\ell}\frac{2\sin\!\big(W'((m-n)\epsilon-\ell P)\big)}{(m-n)\epsilon-\ell P} .\tag{6}$$

Now let $N$ grow. Since $\epsilon\Delta N=2\pi$, two of the three can be chosen.

| held fixed | what goes to zero | the places become | the labels become | the images |
|---|---|---|---|---|
| the length $P$ | $\epsilon$ | a circle of length $P$ | all whole numbers $k$, at $\omega=2\pi k/P$ | stay, one length apart |
| the step $\epsilon$ | $\Delta$ | all whole numbers $n$ | a circle of length $2\pi/\epsilon$ | move off to infinity |
| neither | $\epsilon$ and $\Delta$ | a line | a line | move off to infinity |

34. first row: $W'\to\infty$, each kernel in (6) closes onto a point, and $f$ is a pattern on the circle,

$$f(t)=\sum_{k=-\infty}^{\infty}c_k\,e^{2\pi ikt/P},\qquad c_k=\frac1P\int_0^Pf(t)\,e^{-2\pi ikt/P}\,dt .$$

The group is the turns of a plane by any angle. Its kinds of block are the whole numbers.

35. second row: the group is the whole numbers under addition, steps along a line. Its kinds of block fill a circle: $f(n\epsilon)=\frac1{2\pi}\int_{-\pi/\epsilon}^{\pi/\epsilon}\big[\epsilon\sum_{n'}f(n'\epsilon)\,e^{-i\omega n'\epsilon}\big]e^{i\omega n\epsilon}\,d\omega$.

36. third row: only $\ell=0$ is left and its kernel closes onto a point,

$$f(t)=\frac1{2\pi}\int_{-\infty}^{\infty}\Big[\int_{-\infty}^{\infty}f(t')\,e^{-i\omega t'}\,dt'\Big]e^{i\omega t}\,d\omega ,$$

for functions for which the inner integral exists, for instance those with $\int|f|\,dt$ finite. The circle and the whole numbers are each other's labels; the line is its own.

37. what became of the count. $N=P/\epsilon=2W'/\Delta$: the number of things is the length of the ring over the step, and also the width of the band over the step in label. In every limit $N$ is infinite, but the number of labels below a cut-off $W$ is $2W/\Delta=(W/\pi)\,P$: the count of kinds up to a cut-off, in proportion to the size of the group.

**(E) The class sum of the smallest turns.** Now the rotations of space. In a block of size $d$ a turn by the angle $\theta$ about the direction $\mathbf n$ is $e^{\theta\,(n_xJ_x+n_yJ_y+n_zJ_z)}$, with three generators that obey $J_xJ_y-J_yJ_x=J_z$ and its cyclic turns. The class of that turn is the same angle about every direction. Add the class up, that is, average over $\mathbf n$. Every direction comes with its opposite, so the average of $n_a$ is zero, and the average of $n_an_b$ is $\tfrac13$ for $a=b$ and zero otherwise:

$$\Big\langle e^{\theta\,\mathbf n\cdot\mathbf J}\Big\rangle_{\mathbf n}=1+\frac{\theta^2}{6}\,\big(J_x^2+J_y^2+J_z^2\big)+\dots\tag{7}$$

38. a class sum commutes with everything, at every angle, so the first thing that is left for a small angle does too: the sum of the squares of the generators. It is the class sum of the smallest turns, and like every class sum it is a number on each block:

$$J_x^2+J_y^2+J_z^2=-\frac{d^2-1}{4}\qquad\text{on the block of size }d .\tag{8}$$

39. so "sort by the values of the class sums" (step 3 of section 9) becomes "sort by the value of this one sum": $0,\ -\tfrac34,\ -2,\ -\tfrac{15}4,\dots$ for sizes $1,2,3,4,\dots$, one kind of block for each size. At a finite angle the class is still a number on each block, and by 10.4 the number is the trace over the size. In the basis where the turn about one direction is diagonal its entries are $e^{im\theta}$ for $m=-\tfrac{d-1}2,\dots,\tfrac{d-1}2$, the geometric series of (3) again:

$$\chi_d(\theta)=\frac{\sin(d\,\theta/2)}{\sin(\theta/2)},\qquad\frac{\chi_d(\theta)}d=1-\frac{d^2-1}{24}\,\theta^2+\dots\tag{9}$$

40. the second expression is (7) with (8) in it.

**(F) The rotations by the same four steps.** The classes are the angles $\theta$ from $0$ to $2\pi$, where the turn by $2\pi$ is $-1$ (section 7). The share of the group in the class $\theta$ is the weight $\frac1\pi\sin^2(\theta/2)\,d\theta$; this weight is taken from the standard theory.

*The exact identity.* The rows of the table of traces are perpendicular (10.6), now with the weight:

$$\frac1\pi\int_0^{2\pi}\sin^2\!\frac\theta2\ \chi_d(\theta)\,\chi_{d'}(\theta)\,d\theta=\frac1\pi\int_0^{2\pi}\sin\frac{d\,\theta}2\,\sin\frac{d'\theta}2\,d\theta=\begin{cases}1,&d=d'\\0,&d\neq d' .\end{cases}\tag{10}$$

*The sum over the kinds, with a cut-off.* In place of 10.3, the sum over the blocks up to size $n$:

$$K_n(\theta)=\sum_{d=1}^{n}d\,\chi_d(\theta)=\frac1{\sin(\theta/2)}\sum_{d=1}^{n}d\,\sin\frac{d\,\theta}2$$

$$=-\frac{2}{\sin(\theta/2)}\ \frac{d}{d\theta}\sum_{d=1}^{n}\cos\frac{d\,\theta}2 .\tag{11}$$

41. the last sum is the ring's sum (3) once more. The completeness sum of the rotations is its derivative.

*The limit.* For a function $F$ of the class,

$$\frac1\pi\int_0^{2\pi}\sin^2\!\frac\theta2\ K_n(\theta)\,F(\theta)\,d\theta\ \longrightarrow\ F(0)\qquad(n\to\infty) .\tag{12}$$

With $F(\theta)=e^{\cos(\theta/2)}\,(1+0.3\cos\theta)$, whose value at the identity is $3.53376638$:

| blocks up to size $n$ | the cut-off sum, integrated against $F$ | $\sum_{d\le n}d^2$ |
|---|---|---|
| $4$ | $3.39571311$ | $30$ |
| $8$ | $3.53373263$ | $204$ |
| $16$ | $3.53376638$ | $1496$ |
| $32$ | $3.53376638$ | $11440$ |


*The count.* At the identity $\chi_d=d$, so

$$K_n(0)=\sum_{d=1}^{n}d^2=\frac{n(n+1)(2n+1)}6\ \approx\ \frac{n^3}3 .\tag{13}$$

42. for a finite group $\sum d^2$ is the number of elements. Here it is the number of directions that the blocks up to size $n$ can tell apart, and it grows as the third power of the cut-off. Three is the number of angles it takes to name a rotation. For the turns of a plane the same count is $2K+1$ for labels up to $K$: the first power, for one angle.

**(G) What carries over.**

| finite group | turns by any angle |
|---|---|
| the sum over the elements, divided by $N$ | the integral over the group, with its weight |
| class sums, which commute with everything | class averages; nearest the identity, the sum of the squares of the generators |
| the values of the class sums tell the kinds apart | the value of that one sum does: $-(d^2-1)/4$ |
| a block of size $d$ appears $d$ times | the same: $d^2$ functions for each kind, the entries of its block |
| as many kinds as classes | the classes are a continuum and the kinds a list; the traces are a complete set of functions of the class |
| $\sum d^2=N$ | $\sum_{d\le n}d^2$, a count up to a cut-off, growing as the cut-off to the number of angles |

43. this is for turns. Steps along a line are not turns: their labels are continuous, as in the second and third rows of (D).

**(H) The link to the walk.** The walk's four counters are the block of size two: with the generators $\tfrac12R_{yz}$, $\tfrac12R_{zx}$, $\tfrac12R_{xy}$ the sum of the squares is $-\tfrac34$, the value of (8) for $d=2$ (the three obey the rule of (E) up to the sign of one of them). The walker's own turn is a turn in a plane, so its kinds of block are the whole numbers of (D), first row: how many times the walker's hand turns for one turn of the field's angle. The walk uses the kind $1$, and $-1$ for the added four counted with the other sign.

A warning on names: the sum of squares in (8) is called the Casimir operator. It has nothing to do with the Casimir energy except the person. What does connect to the Casimir energy is (5): a sum is its integral plus the terms that go round and come back. Part (I) carries that out.

**(I) The partition function, the Casimir energy and zeta.** Units $\hbar=c=1$ in this part.

*The picture first.* Give every label a weight that falls off with its frequency. The sum of the weights is the partition function of the ring. The rule "a sum over labels is an integral plus its images" writes the same function as a sum over windings: the number of times something goes round the ring before it comes back. The energy of all the labels is an integral of the partition function. In it the term with no winding grows with the length of the ring, and the terms that wind are what is left: the Casimir energy. Seen from the side of the labels the same numbers are values of zeta.

*The rule in general.* (5) and item 33, for any weight $g(\omega)$ on the labels $\omega_k=2\pi k/P$:

$$\sum_{k=-\infty}^{\infty}g\Big(\frac{2\pi k}{P}\Big)=\frac{P}{2\pi}\sum_{\ell=-\infty}^{\infty}\int_{-\infty}^{\infty}g(\omega)\,e^{i\omega\,\ell P}\,d\omega .\tag{14}$$

44. the term $\ell=0$ is the integral over the labels. The term $\ell$ is the same integral with a wave that fits $\ell$ lengths of the ring.

*The partition function of the ring.* With the weight $g=e^{-\tau\omega^2}$, and $\int e^{-\tau\omega^2+i\omega x}\,d\omega=\sqrt{\pi/\tau}\;e^{-x^2/4\tau}$:

$$Z_P(\tau)=\sum_{k}e^{-\tau(2\pi k/P)^2}=\frac{P}{\sqrt{4\pi\tau}}\sum_{\ell}e^{-\ell^2P^2/4\tau} .\tag{15}$$

45. two forms of one function: a sum over the labels, and a sum over the windings. On the right, $\frac1{\sqrt{4\pi\tau}}e^{-x^2/4\tau}$ is Part I's spreading after a time $\tau$: the chance that a diffusing point is found a distance $x$ away. So $Z_P(\tau)/P$ is the chance that it is back where it started, added over the number of times $\ell$ it has gone round. For large $\tau$ the left form is short (only $k=0$ counts); for small $\tau$ the right form is (only $\ell=0$).

*With the flip rate.* The energy of label $k$ is $E_k=\sqrt{a^2+\omega_k^2}$, and the weight $e^{-\tau E_k^2}$ only gains a factor:

$$\sum_ke^{-\tau E_k^2}=e^{-\tau a^2}\,Z_P(\tau) .\tag{16}$$

*An energy as an integral of its weight.* For any $E>0$,

$$E=-\frac1{2\sqrt\pi}\int_0^\infty\tau^{-3/2}\,\big(e^{-\tau E^2}-1\big)\,d\tau .\tag{17}$$

46. the $-1$ does not depend on $E$. It is left out below: it adds the same constant to every label, and in the form of windings that falls in the term $\ell=0$ only.

*The energy of all the labels.* Put (17) on each label, use (16), and write $Z_P$ in its form of windings. The term $\ell$ is

$$-\frac1{2\sqrt\pi}\int_0^\infty\tau^{-3/2}\,e^{-\tau a^2}\,\frac{P}{\sqrt{4\pi\tau}}\,e^{-\ell^2P^2/4\tau}\,d\tau=-\frac{P}{4\pi}\int_0^\infty\frac{d\tau}{\tau^{2}}\ e^{-a^2\tau-\ell^2P^2/4\tau}$$

$$=-\frac{a}{\pi\,|\ell|}\,K_1\big(a\,|\ell|\,P\big)\qquad(\ell\neq0),$$

so that

$$\sum_k\sqrt{a^2+\omega_k^2}\ =\ \underbrace{\frac{P}{2\pi}\int\sqrt{a^2+\omega^2}\;d\omega}_{\ell=0:\ \text{grows with }P}\ \ \underbrace{-\ \frac{2a}{\pi}\sum_{\ell=1}^{\infty}\frac{K_1(a\,\ell\,P)}{\ell}}_{\text{the windings}} .\tag{18}$$

47. the first term is infinite and proportional to $P$: the energy of the line, per length, times the length. The second is finite. Checked with a smooth cut-off: for $a=1$, $P=3$ the sum over the labels minus its integral is $-0.02600$, and the windings give $-0.02600$.

48. this is the Casimir document's formula. Its integral $\int\frac{d\tau}{\tau}\,\frac{e^{-(2nL)^2/4\tau-m^2\tau}}{4\pi\tau}$, for a loop that winds $n$ times between walls a distance $L$ apart, is the integral above with $P=2L$, $\ell=n$, $a=m$. Its energy $E(L)=-\frac m{2\pi}\sum_n\frac{K_1(2nmL)}n$ is a quarter of the windings of (18): a half because each label carries half its frequency, and a half because the walls keep half the labels of the ring of length $2L$.

49. the flip rate in the windings: $K_1(a\ell P)$ falls off as $e^{-a\ell P}$. A loop that goes round a ring longer than $1/a$ is suppressed. With no flip rate, $K_1(z)\to1/z$ and the windings are $-\frac{2}{\pi P}\sum_\ell\frac1{\ell^2}=-\frac{\pi}{3P}$.

*Zeta.* With no flip rate the labels can be summed directly. For the labels other than $k=0$,

$$\sum_{k\neq0}|\omega_k|^{-s}=2\Big(\frac{P}{2\pi}\Big)^{s}\zeta(s),\tag{19}$$

$$\text{at }s=-1:\qquad2\cdot\frac{2\pi}{P}\cdot\zeta(-1)=-\frac{\pi}{3P}\qquad\text{with }\zeta(-1)=-\tfrac1{12} .$$

50. the same number as item 49, found on the other side. There it came from $\sum1/\ell^2=\zeta(2)=\pi^2/6$ over the windings; here from $\zeta(-1)$ over the labels. The two are tied, $\zeta(-1)=-\zeta(2)/2\pi^2$, and what ties them is (15): zeta is an integral of the partition function,

$$\sum_{k\neq0}|\omega_k|^{-s}=\frac1{\Gamma(s/2)}\int_0^\infty\tau^{s/2-1}\,\big(Z_P(\tau)-1\big)\,d\tau ,\tag{20}$$

for $s>1$ and elsewhere by continuation, and exchanging the two forms of $Z_P$ exchanges $s$ with $1-s$.

*One label.* Make the ring short, $P\ll1/a$. In (15) every label but $k=0$ is then far above $a$, and for any fixed $\tau$ the partition function with the flip rate is $e^{-\tau a^2}$: one label, of energy $a$. That is the walker at rest, with nothing but its clock. The flip rate is not produced by this limit; it is what the limit leaves.

*The rotations.* The same function for the rotations: each kind of block with its count $d^2$ and the weight of its value (8). By the same rule,

$$Z_{\rm rot}(\tau)=\sum_{d=1}^{\infty}d^2\,e^{-\tau(d^2-1)/4}=e^{\tau/4}\ \frac{2\sqrt\pi}{\tau^{3/2}}\sum_{\ell=-\infty}^{\infty}\Big(1-\frac{8\pi^2\ell^2}{\tau}\Big)\,e^{-4\pi^2\ell^2/\tau} .\tag{21}$$

51. sizes on the left, windings on the right: a winding is a turn by $4\pi\ell$, the angle after which a block of every size is back. For small $\tau$ the term $\ell=0$ is $2\sqrt\pi/\tau^{3/2}$: the power is half the number of angles, the same statement as the count (13). What (21) is for the walk is not established: the labels the walk sums over are those of positions, as in (18), not the sizes of the blocks of the rotations.

*A remark on counting.* The geometric series behind (3) and (9), $q^{-(d-1)/2}+\dots+q^{(d-1)/2}$, is called a $q$-number. At $q=e^{i\theta}$ it is the trace of the block of size $d$, and at $q\to1$, the identity, it is the plain count $d$. The step that was replaced in (C), $1-e^{\Delta D}\to-\Delta D$, is this limit taken inside the sum. "The field with one element" is a name for doing everything at $q=1$, where matrices become permutations; it is a programme, not a settled object, and nothing here depends on it.

## 13. The same things in the standard names

| in this note | standard name |
|---|---|
| a group acting on its own elements as permutations | Cayley's theorem; the regular representation |
| a permutation is made of swaps; a length-keeping matrix is made of mirrors | transpositions generate the symmetric group; the Cartan–Dieudonné theorem |
| one more line carrying $\det g$ | the embedding of $O(n)$ in $SO(n+1)$; the sign (determinant) representation |
| patterns with $k$ waves make every shift diagonal | characters of a commuting group; the discrete Fourier transform |
| a block | an irreducible representation |
| what commutes with everything only mixes identical blocks | Schur's lemma |
| a block of size $d$ appears $d$ times; $\sum d^2$ is the size of the group | decomposition of the regular representation |
| kinds of block = classes; rows of traces perpendicular | the character table and its orthogonality |
| the entries of the block matrices as the new basis | the Fourier transform on a group (Peter–Weyl) |
| the group of mirrors is twice the group of turns; the full turn is $-1$ | the Pin and Spin groups; the double cover |
| steps 1 to 7 | decomposition of the regular representation into isotypic components and irreducible ones; related to the Burnside–Dixon method for the table of traces |
| a block that cannot be split is filled completely | Burnside's theorem (the double commutant) |
| the number $\lambda_{ij}$ of a class sum on a kind | the central character |
| 10.5, 10.6 | the projection formula; orthogonality of characters |
| 10.7 | the Fourier transform on a finite group, its inverse, the Plancherel formula, the convolution theorem |
| a sum is its integral plus its images, (5) and item 33 | the Poisson summation formula; aliasing |
| the three limits of (D) | Fourier series, the discrete-time Fourier transform, the Fourier integral; Pontryagin duality |
| the sum of the squares of the generators, (8) | the Casimir operator of the rotations; $j(j+1)$ with $d=2j+1$ |
| the traces (9) and the weight of (F) | the characters of $SU(2)$; the Haar measure on classes (Weyl's integration formula) |
| (11), (12) | completeness of the characters (Peter–Weyl) |
| the count (13) | Weyl's law |
| (14) | the Poisson summation formula |
| (15) and its two forms | the heat kernel on a circle; the Jacobi theta function and its inversion; the method of images |
| (17), (18) | proper-time (heat-kernel) regularization of the vacuum energy; the Casimir energy of a massive field on a circle |
| (19), (20) | zeta-function regularization; the Mellin transform; the functional equation of the Riemann zeta function |
| (21) | the heat kernel on $SU(2)$, as a sum over characters and as a sum over geodesics |
| $q$-numbers | $q$-analogues; quantum dimensions |

(`scripts/groups_as_permutations.py`: 1 for sections 1–2; 2 for section 3; 3 for section 4; 4–6 for section 5; 7–8 for section 6; 9 for section 7. `scripts/group_blocks_algorithm.py`: sections 9–11; each line of its output is one group, with the largest error over 10.1–10.7. `scripts/fixed_angle_to_any_angle.py`: 1 for (1), (2); 2 for (3); 3 for (4), (5); 4 for item 33; 5 for (7)–(9); 6 for (10); 7 for (12), (13); 8 for (H); 9 for (15); 10 for (17); 11 for (18) and item 48; 12 for (19), (20); 13 for (21); 14 for the remark on counting.)
