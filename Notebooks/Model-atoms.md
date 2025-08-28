```wl
In[]:= \[AliasDelimiter]
```

```wl
In[]:= 
```

# Solving The Schrödinger Equation for Finite Potential: from atoms to molecules to solids

Authors: Angel Salazar 
School of Physical Sciences and Nanotechnology

### Single Finite Potential Well: model of a single atom

We want to solve the time independent Schrödinger equation:
$\frac{-\hbar ^2}{2m}$ $\frac{\partial ^2\psi (x)}{\partial x}$ V(x)ψ(x) = E ψ(x)

Let's define the finite potential.  use +pw+:

![0gvuglbmam2iv](img/0gvuglbmam2iv.png)

![1na3obejdemjj](img/1na3obejdemjj.png)

```wl
In[]:= Plot[v[x] /. {v0 -> vpot, a -> la}, {x, -la - 0.3, la + 0.3}, 
   AxesStyle -> {"x", "V(x)"}, 
   PlotStyle -> {Red, Thick}, 
   PlotRange -> {Automatic, {-0.2, vpot + 1.5}}, 
   AxesStyle -> Arrowheads[{0.0, 0.03}], 
   Filling -> Axis 
  ]
```

![0ntpxb1zh1ddo](img/0ntpxb1zh1ddo.png)

#### Solving the Schrödinger equation in a . u . 

In atomic units, ℏ -> 1, m -> 1, e -> 1; then the lengths are given in Bohr, 1bohr= 0.53 Å, and the energy in Hartree $E_h$ = 27.21 eV. 
Let's define α = $\sqrt{2(\text{v0} - \mathcal{E})}$, β = $\sqrt{2 \mathcal{E}}$,   (ℰ is +scE+)

![08az2jueo4sbh](img/08az2jueo4sbh.png)

![1f67s4tz6sbir](img/1f67s4tz6sbir.png)

![1qrjt5ner8k04](img/1qrjt5ner8k04.png)

![1ql7zaha7bnlq](img/1ql7zaha7bnlq.png)

![0zaxipanmaqgm](img/0zaxipanmaqgm.png)

![0zi533y8qmjto](img/0zi533y8qmjto.png)

![10b261ngxgdiw](img/10b261ngxgdiw.png)

![0jdeepjf7l23u](img/0jdeepjf7l23u.png)

![0gjksx8xanmq3](img/0gjksx8xanmq3.png)

Let's define the wavefunction  φ, φ is +j+:

![07nyzvkcmbljy](img/07nyzvkcmbljy.png)

```wl
In[]:= 
```

Now we need to consider the boundary conditions

![0zes6cyw5fqlx](img/0zes6cyw5fqlx.png)

![1oall2lhm5py2](img/1oall2lhm5py2.png)

![0uxxvogu1q5dj](img/0uxxvogu1q5dj.png)

![0p7s1hyn8x11a](img/0p7s1hyn8x11a.png)

![13vswao35eg61](img/13vswao35eg61.png)

![0s29jcfqjmmgh](img/0s29jcfqjmmgh.png)

![0vawdlda6p8x5](img/0vawdlda6p8x5.png)

```wl
In[]:= mat // MatrixForm
```

|  |  |  |  |
| - | - | - | - |
| -Power- | -Times- | -Sin- | 0 |
| -Times- | -Times- | -Times- | 0 |
| 0 | -Cos- | -Sin- | -Times- |
| 0 | -Times- | -Times- | -Times- |

![0g32t1sh57fj5](img/0g32t1sh57fj5.png)

![15iokcbn8btyk](img/15iokcbn8btyk.png)

```wl
In[]:= Plot[Det[mat2], {\[ScriptCapitalE], 0, vpot}]
```

![05at05k2t7xpg](img/05at05k2t7xpg.png)

```wl
In[]:= fig = Plot[Det[mat2], {\[ScriptCapitalE], 0, vpot}, PlotRange -> {-10^-20, 10^-20}]
```

![085gtgosabwo0](img/085gtgosabwo0.png)

```wl
In[]:= Short[fig, 5]
```

![0y1qo6sw6jcos](img/0y1qo6sw6jcos.png)

```wl
In[]:= fig[[1, 1, 1, 1, 1, 3, 1]][[2 ;; -1]]
```

```wl
Out[]= {Line[{{0.748476, 1.*10^-20}, {0.748476, -1.*10^-20}}], Line[{{2.9027, -1.*10^-20}, {2.9027, 1.*10^-20}}], Line[{{5.92702, 1.*10^-20}, {5.92702, -1.*10^-20}}]}
```

```wl
In[]:= approxroot1 = fig[[1, 1, 1, 1, 1, 3, 1]][[2 ;; -1]] // Map[#[[1, 1, 1]] &, #] & // Sort
```

```wl
Out[]= {0.748476, 2.9027, 5.92702}
```

```wl
In[]:= rval1 = approxroot1 // Map[FindRoot[Det[mat2], {\[ScriptCapitalE], #}] &, #] & // Map[\[ScriptCapitalE] /. # &, #] &
```

```wl
Out[]= {0.7495, 2.9032, 5.92714}
```

```wl
In[]:= fig[[1, 1, 1]]
```

![0wo9x42wz8avp](img/0wo9x42wz8avp.png)

```wl
In[]:= g = Normal[fig];  
   
  (*extract the x of every vertical line the plot drew (your zero crossings)*)
 aproxroot = Sort@Cases[g, Line[{{x_, y1_}, {x_, y2_}}] :> x, Infinity]
 
```

```wl
Out[]= {0.748476, 2.9027, 5.92702}
```

```wl
In[]:= rval = aproxroot // Map[FindRoot[Det[mat2], {\[ScriptCapitalE], #}] &, #] & // Map[\[ScriptCapitalE] /. # &, #] &
```

```wl
Out[]= {0.7495, 2.9032, 5.92714}
```

#### Level 1

![0o1u4t6o5vknb](img/0o1u4t6o5vknb.png)

![0tpqkqpy35etz](img/0tpqkqpy35etz.png)

![1l7ls47ko3uso](img/1l7ls47ko3uso.png)

![0rl6c6q4e2it8](img/0rl6c6q4e2it8.png)

![1m1y9hrm5hzuf](img/1m1y9hrm5hzuf.png)

```wl
Out[]= {0.705376, 0.0699299, 0, 0.705376}
```

We need to find the normalization constant for φ(x) , ∞

![0rv40uoab70bu](img/0rv40uoab70bu.png)

```wl
Out[]= 0.00633217
```

![02h0fnbs2oulp](img/02h0fnbs2oulp.png)

![0v6zefkibmond](img/0v6zefkibmond.png)

```wl
In[]:= 
```

![104zpzzvn20hj](img/104zpzzvn20hj.png)

```wl
Out[]= 1.
```

What is the probability of finding one electron when $x>\frac{a}{2}$?

![0rnjd5nrpiunc](img/0rnjd5nrpiunc.png)

```wl
Out[]= 0.0131291
```

Let' s find the Δx, then we need the expected value of 〈x〉 and $\left\langle x^2\right\rangle$

![0r3nezpnb8wvh](img/0r3nezpnb8wvh.png)

![0fof3kb6xn718](img/0fof3kb6xn718.png)

```wl
Out[]= 0
```

![15rlrasqlou4f](img/15rlrasqlou4f.png)

```wl
Out[]= 0.228641
```

Then the Δx is :

![06utkzqrrftit](img/06utkzqrrftit.png)

```wl
Out[]= 0.478165
```

Now let' s find the uncertainty in the momentum Δp :

![11yn8x0lzw1p0](img/11yn8x0lzw1p0.png)

![0l7gki74qejxd](img/0l7gki74qejxd.png)

![00w504jdl0vuq](img/00w504jdl0vuq.png)

```wl
Out[]= 0
```

![1by4rdnewcbdt](img/1by4rdnewcbdt.png)

```wl
Out[]= 1.15764
```

Then Δp is 

![1721j7bu648va](img/1721j7bu648va.png)

```wl
Out[]= 1.07594
```

Therefore the Heisengerg uncertainty Δx Δp is:

```wl
In[]:= \[CapitalDelta]x \[CapitalDelta]p
```

```wl
Out[]= 0.514476
```

What is the expected value of the kinetic energy k = $\frac{-\hbar ^2}{2m}$

![0i6my177679yn](img/0i6my177679yn.png)

```wl
Out[]= 0.578822
```

What is the expected potential energy for this level?

![0crs0qx6dgjmi](img/0crs0qx6dgjmi.png)

```wl
Out[]= 0.170678
```

Now, check this out: The total energy of the level

![1p9qrgrbn726z](img/1p9qrgrbn726z.png)

```wl
Out[]= 0.7495
```

```wl
In[]:= \[ScriptCapitalE] /. \[ScriptCapitalE] -> rval[[nval]]
```

```wl
Out[]= 0.7495
```

#### Level 2

![10cjh3eh4ylph](img/10cjh3eh4ylph.png)

```wl
Out[]= {-0.705261, 0, 0.0722025, 0.705261}
```

We need to find the normalization constant for φ(x) , ∞

![0d32n2y29lj3k](img/0d32n2y29lj3k.png)

```wl
Out[]= 0.00715691
```

![1dajukkqjdqrd](img/1dajukkqjdqrd.png)

![1erpo38mfam0o](img/1erpo38mfam0o.png)

```wl
In[]:= 
```

![1b2fmf09yr5ns](img/1b2fmf09yr5ns.png)

```wl
Out[]= 1.
```

What is the probability of finding one electron when $x>\frac{a}{2}$?

![1fkyfam0j6v8r](img/1fkyfam0j6v8r.png)

```wl
Out[]= 0.0606513
```

Let' s find the Δx, then we need the expected value of 〈x〉 and $\left\langle x^2\right\rangle$

![1hqkrb5bkd2d9](img/1hqkrb5bkd2d9.png)

![01rqlpsgbj5i0](img/01rqlpsgbj5i0.png)

```wl
Out[]= 0
```

![0kt51d173miwz](img/0kt51d173miwz.png)

```wl
Out[]= 0.548414
```

Then the Δx is :

![0hbddmstiadfo](img/0hbddmstiadfo.png)

```wl
Out[]= 0.74055
```

Now let' s find the uncertainty in the momentum Δp :

![1ul32162v1f8n](img/1ul32162v1f8n.png)

![141v033ci7i52](img/141v033ci7i52.png)

```wl
Out[]= 0
```

![101vv2tnlidgx](img/101vv2tnlidgx.png)

```wl
Out[]= 4.22948
```

Then Δp is 

![1l5glr7pfmtr0](img/1l5glr7pfmtr0.png)

```wl
Out[]= 2.05657
```

Therefore the Heisengerg uncertainty Δx Δp is:

```wl
In[]:= \[CapitalDelta]x \[CapitalDelta]p
```

```wl
Out[]= 1.52299
```

What is the expected value of the kinetic energy k = $\frac{-\hbar ^2}{2m}$

![1q895gxw1u8db](img/1q895gxw1u8db.png)

```wl
Out[]= 2.11474
```

What is the expected potential energy for this level?

![1u4kaiwtqzfv2](img/1u4kaiwtqzfv2.png)

```wl
Out[]= 0.788467
```

Now, check this out: The total energy of the level

![0qsm21vy5hxn2](img/0qsm21vy5hxn2.png)

```wl
Out[]= 2.9032
```

```wl
In[]:= \[ScriptCapitalE] /. \[ScriptCapitalE] -> rval[[nval]]
```

```wl
Out[]= 2.9032
```

#### Level 3

![0jgq305if6r0j](img/0jgq305if6r0j.png)

```wl
Out[]= {0.685361, -0.24609, 0, 0.685361}
```

We need to find the normalization constant for φ(x) , ∞

![0jnmr4hwgs5ib](img/0jnmr4hwgs5ib.png)

```wl
Out[]= 0.117138
```

![07xr1os4qfa41](img/07xr1os4qfa41.png)

![0112riy05l9yw](img/0112riy05l9yw.png)

```wl
In[]:= 
```

![1olsybhsut1uf](img/1olsybhsut1uf.png)

```wl
Out[]= 1.
```

What is the probability of finding one electron when $x>\frac{a}{2}$?

![0vrgn2wchlqrf](img/0vrgn2wchlqrf.png)

```wl
Out[]= 0.220218
```

Let' s find the Δx, then we need the expected value of 〈x〉 and $\left\langle x^2\right\rangle$

![1bu89w6j5qilc](img/1bu89w6j5qilc.png)

![0avixi13td5yx](img/0avixi13td5yx.png)

![0bjtvr4dz72uj](img/0bjtvr4dz72uj.png)

```wl
Out[]= 0
```

![1uv28ckop0fgb](img/1uv28ckop0fgb.png)

```wl
Out[]= 1.27519
```

Then the Δx is :

![1bea7mzvahb8n](img/1bea7mzvahb8n.png)

```wl
Out[]= 1.12924
```

Now let' s find the uncertainty in the momentum Δp :

![1rgbcmaalu2h3](img/1rgbcmaalu2h3.png)

![1dkhz2nuiforn](img/1dkhz2nuiforn.png)

```wl
Out[]= 0
```

![0it2uofboihh8](img/0it2uofboihh8.png)

```wl
Out[]= 6.12863
```

Then Δp is 

![0zw2w79d2p33q](img/0zw2w79d2p33q.png)

```wl
Out[]= 2.47561
```

Therefore the Heisengerg uncertainty Δx Δp is:

```wl
In[]:= \[CapitalDelta]x \[CapitalDelta]p
```

```wl
Out[]= 2.79556
```

What is the expected value of the kinetic energy k = $\frac{-\hbar ^2}{2m}$

![0v4idb2bmm7p6](img/0v4idb2bmm7p6.png)

```wl
Out[]= 3.06431
```

What is the expected potential energy for this level?

![0jy2jg6g7nnfa](img/0jy2jg6g7nnfa.png)

```wl
Out[]= 2.86283
```

Now, check this out: The total energy of the level

![1a71no3scn7nx](img/1a71no3scn7nx.png)

```wl
Out[]= 5.92714
```

```wl
In[]:= \[ScriptCapitalE] /. \[ScriptCapitalE] -> rval[[nval]]
```

```wl
Out[]= 5.92714
```

#### Summary of the plots

```wl
In[]:= Show @@ Table[figure[i], {i, 3}]
```

![0m073vi80btag](img/0m073vi80btag.png)

What is the total energy of the system if it has 3 E?

```wl
In[]:= 
```

![1s0k9s4kxum5e](img/1s0k9s4kxum5e.png)

```wl
Out[]= 4.40221
```

![13rct37euxpez](img/13rct37euxpez.png)

```wl
Out[]= 3.27238
```

![13q5ehe6u6ycr](img/13q5ehe6u6ycr.png)

```wl
Out[]= 1.12982
```

```wl
In[]:= 
```

### Two Finite Potential Well: model of a molecule

We want to solve the time independent Schrödinger equation:
$\frac{-\hbar ^2}{2m}$ $\frac{\partial ^2\psi (x)}{\partial x}$ V(x)ψ(x) = E ψ(x)

#### Define the potential

Let's define the finite potential.  use +pw+:

![1o0zr1zp8dq6a](img/1o0zr1zp8dq6a.png)

![1ulukt1okmsgc](img/1ulukt1okmsgc.png)

```wl
In[]:= Plot[v[x] /. {v0 -> vpot, a -> la, b -> lb}, {x, -la - lb - 0.3, la + la + 0.3}, 
   AxesStyle -> {"x", "V(x)"}, 
   PlotStyle -> {Red, Thick}, 
   PlotRange -> {Automatic, {-0.2, vpot + 1.5}}, 
   AxesStyle -> Arrowheads[{0.0, 0.03}], 
   Filling -> Axis 
  ]
```

![1o3r6zhkflcll](img/1o3r6zhkflcll.png)

#### Solving the Schrödinger equation in a . u . 

In atomic units, ℏ -> 1, m -> 1, e -> 1; then the lengths are given in Bohr, 1bohr= 0.53 Å, and the energy in Hartree $E_h$ = 27.21 eV. 
Let's define α = $\sqrt{2(\text{v0} - \mathcal{E})}$, β = $\sqrt{2 \mathcal{E}}$,   (ℰ is +scE+)

![1xaagxbrcu78e](img/1xaagxbrcu78e.png)

![02pxihbis7vm0](img/02pxihbis7vm0.png)

![1ja2wdsg28s3c](img/1ja2wdsg28s3c.png)

![1dheq3sh4fy0h](img/1dheq3sh4fy0h.png)

![04021i41mw00h](img/04021i41mw00h.png)

![0i93n5p6776ks](img/0i93n5p6776ks.png)

![1bdolta09bbhp](img/1bdolta09bbhp.png)

![0mcr22zuxbztk](img/0mcr22zuxbztk.png)

![0x5f6ncjwe4fj](img/0x5f6ncjwe4fj.png)

Let's define the wavefunction  φ, φ is +j+:

![07szzh60ies2r](img/07szzh60ies2r.png)

```wl
In[]:= 
```

Now we need to consider the boundary conditions

![1ce18krlkaly3](img/1ce18krlkaly3.png)

![0qpprod1qplez](img/0qpprod1qplez.png)

![19t20clpqcb9r](img/19t20clpqcb9r.png)

![0qwik0txi2sz3](img/0qwik0txi2sz3.png)

![0so86dwhqglxd](img/0so86dwhqglxd.png)

![1u0j8ku951erl](img/1u0j8ku951erl.png)

![082z56qyaan2d](img/082z56qyaan2d.png)

![1fhrsjolv8lif](img/1fhrsjolv8lif.png)

![18bd6rmn5paql](img/18bd6rmn5paql.png)

![0oc3wk8g848g1](img/0oc3wk8g848g1.png)

![1evu5627vz8nz](img/1evu5627vz8nz.png)

```wl
In[]:= mat // MatrixForm
```

|  |  |  |  |  |  |  |  |
| - | - | - | - | - | - | - | - |
| -Power- | -Times- | -Times- | 0 | 0 | 0 | 0 | 0 |
| -Times- | -Times- | -Times- | 0 | 0 | 0 | 0 | 0 |
| 0 | -Cos- | -Times- | -Times- | -Times- | 0 | 0 | 0 |
| 0 | -Times- | -Times- | -Times- | -Times- | 0 | 0 | 0 |
| 0 | 0 | 0 | -Power- | -Power- | -Times- | -Times- | 0 |
| 0 | 0 | 0 | -Times- | -Times- | -Times- | -Times- | 0 |
| 0 | 0 | 0 | 0 | 0 | -Cos- | -Sin- | -Times- |
| 0 | 0 | 0 | 0 | 0 | -Times- | -Times- | -Times- |

![1wo18qree3ndg](img/1wo18qree3ndg.png)

![1spk4q24plurd](img/1spk4q24plurd.png)

```wl
In[]:= Plot[Det[mat2], {\[ScriptCapitalE], 0, vpot}]
```

![011a3r36bzor9](img/011a3r36bzor9.png)

```wl
In[]:= fig = Plot[Det[mat2], {\[ScriptCapitalE], 0, vpot}, PlotRange -> {-10^-20, 10^-20}]
```

![1m59dzi9ahk4o](img/1m59dzi9ahk4o.png)

```wl
In[]:= fig[[1, 1, 1, 1, 1, 3, 1]][[2 ;; -1]]
```

```wl
Out[]= {Line[{{5.78811, 1.*10^-20}, {5.78811, -1.*10^-20}}], Line[{{0.739199, 1.*10^-20}, {0.739199, -1.*10^-20}}], Line[{{6.19469, -1.*10^-20}, {6.19469, 1.*10^-20}}], Line[{{0.7595, -1.*10^-20}, {0.7595, 1.*10^-20}}], Line[{{2.84476, 1.*10^-20}, {2.84476, -1.*10^-20}}], Line[{{2.96416, -1.*10^-20}, {2.96416, 1.*10^-20}}]}
```

```wl
In[]:= approxroot = fig[[1, 1, 1, 1, 1, 3, 1]][[2 ;; -1]] // Map[#[[1, 1, 1]] &, #] & // Sort
```

```wl
Out[]= {0.739199, 0.7595, 2.84476, 2.96416, 5.78811, 6.19469}
```

```wl
In[]:= rval = approxroot // Map[FindRoot[Det[mat2], {\[ScriptCapitalE], #}] &, #] & // Map[\[ScriptCapitalE] /. # &, #] &
```

```wl
Out[]= {0.739194, 0.759528, 2.84474, 2.96417, 5.78909, 6.19471}
```

```wl
In[]:= fig[[1, 1, 1]]
```

![05q3a8cpms7j6](img/05q3a8cpms7j6.png)

#### Level 1

![01ex9e1qzl89h](img/01ex9e1qzl89h.png)

![0zr2rc9lwqgt1](img/0zr2rc9lwqgt1.png)

![0kt0uvwh7qots](img/0kt0uvwh7qots.png)

![0y3fv36064gxq](img/0y3fv36064gxq.png)

![07nahf37c9ahq](img/07nahf37c9ahq.png)

```wl
Out[]= {0.707107, -0.000103738, -0.000420087, 0.0000274377, 0.0000274377, -0.000103738, 0.000420087, 0.707107}
```

We need to find the normalization constant for φ(x) , ∞

![15gu7epci2ddt](img/15gu7epci2ddt.png)

```wl
Out[]= 4.8918*10^-7
```

![06o1xyjlfis8j](img/06o1xyjlfis8j.png)

![1qt5aqpqmkun0](img/1qt5aqpqmkun0.png)

**Here we can see the bounding state** 

For this first energy level we can see that there is no nodes, since we have two atoms the energy level splits, it is forming a bounding level. This is the origin of the covalent bond which is related with the tunneling effect, since in the forbidden region is non zero, so the electron can be found with probabilities to be in the left region or right region. The blue line is the wave function. 

```wl
In[]:= 
```

![1quqb57lw14xj](img/1quqb57lw14xj.png)

```wl
Out[]= 1.
```

What is the probability of finding one electron when $x>a+\frac{b}{2}$?

![0qqqsgiap540w](img/0qqqsgiap540w.png)

```wl
Out[]= 0.166753
```

Let' s find the Δx, then we need the expected value of 〈x〉 and $\left\langle x^2\right\rangle$

![0lwq11cpjt6c1](img/0lwq11cpjt6c1.png)

![1v100be0qw483](img/1v100be0qw483.png)

![14shdahs75w7v](img/14shdahs75w7v.png)

```wl
Out[]= 1.00616*10^-10
```

![03gob1zynxqjt](img/03gob1zynxqjt.png)

```wl
Out[]= 2.44961
```

Then the Δx is :

![10jhjymlqzqdh](img/10jhjymlqzqdh.png)

```wl
Out[]= 1.56512
```

Now let' s find the uncertainty in the momentum Δp :

![17exb9yrnxt8b](img/17exb9yrnxt8b.png)

![1437p9ebl86pw](img/1437p9ebl86pw.png)

```wl
Out[]= 0
```

![09i1x1u98wtwn](img/09i1x1u98wtwn.png)

```wl
Out[]= 1.09625
```

Then Δp is 

![18f32ust2gtpq](img/18f32ust2gtpq.png)

```wl
Out[]= 1.04702
```

Therefore the Heisengerg uncertainty Δx Δp is:

```wl
In[]:= \[CapitalDelta]x \[CapitalDelta]p
```

```wl
Out[]= 1.63872
```

What is the expected value of the kinetic energy k = $\frac{-\hbar ^2}{2m}$

![1q6msdxsn4ntx](img/1q6msdxsn4ntx.png)

```wl
Out[]= 0.548126
```

What is the expected potential energy for this level?

![1raee8byxvt64](img/1raee8byxvt64.png)

```wl
Out[]= 0.191068
```

Now, check this out: The total energy of the level

![1txz6ubrnor87](img/1txz6ubrnor87.png)

```wl
Out[]= 0.739194
```

```wl
In[]:= \[ScriptCapitalE] /. \[ScriptCapitalE] -> rval[[nval]]
```

```wl
Out[]= 0.739194
```

#### Level 2

![13bnjw1aoc816](img/13bnjw1aoc816.png)

```wl
Out[]= {-0.707107, 0.000123297, 0.000415411, -0.0000265255, 0.0000265255, -0.000123297, 0.000415411, 0.707107}
```

We need to find the normalization constant for φ(x) , ∞

![0vfh4agb71yxy](img/0vfh4agb71yxy.png)

```wl
Out[]= 4.82117*10^-7
```

![08ty5woct7d9k](img/08ty5woct7d9k.png)

![15de86u99k3li](img/15de86u99k3li.png)

**Here we can see the bounding state** 

For this first energy level we can see that there is no nodes, since we have two atoms the energy level splits, it is forming a bounding level. This is the origin of the covalent bond which is related with the tunneling effect, since in the forbidden region is non zero, so the electron can be found with probabilities to be in the left region or right region. The blue line is the wave function. 

```wl
In[]:= 
```

![0zuw82k389jrz](img/0zuw82k389jrz.png)

```wl
Out[]= 1.
```

What is the probability of finding one electron when $x>a+\frac{b}{2}$?

![0kzwj8e0ir64u](img/0kzwj8e0ir64u.png)

```wl
Out[]= 0.15074
```

Let' s find the Δx, then we need the expected value of 〈x〉 and $\left\langle x^2\right\rangle$

![1s3g63ea0ci0k](img/1s3g63ea0ci0k.png)

![09visbwn27x7n](img/09visbwn27x7n.png)

![1eyc97z4hk8g5](img/1eyc97z4hk8g5.png)

```wl
Out[]= 1.32334*10^-10
```

![1mnma5qlicshw](img/1mnma5qlicshw.png)

```wl
Out[]= 2.50693
```

Then the Δx is :

![19w27c0kyqx2n](img/19w27c0kyqx2n.png)

```wl
Out[]= 1.58333
```

Now let' s find the uncertainty in the momentum Δp :

![19phaecvowbbu](img/19phaecvowbbu.png)

![0e6yftko57coa](img/0e6yftko57coa.png)

```wl
Out[]= 0
```

![0ref2jqh4jcsl](img/0ref2jqh4jcsl.png)

```wl
Out[]= 1.21675
```

Then Δp is 

![0j8dwcp8tqqm2](img/0j8dwcp8tqqm2.png)

```wl
Out[]= 1.10306
```

Therefore the Heisengerg uncertainty Δx Δp is:

```wl
In[]:= \[CapitalDelta]x \[CapitalDelta]p
```

```wl
Out[]= 1.74651
```

What is the expected value of the kinetic energy k = $\frac{-\hbar ^2}{2m}$

![0asdxrksackfu](img/0asdxrksackfu.png)

```wl
Out[]= 0.608376
```

What is the expected potential energy for this level?

![0f1mfxmnhe9o4](img/0f1mfxmnhe9o4.png)

```wl
Out[]= 0.151152
```

Now, check this out: The total energy of the level

![0j9zfsbjor1s6](img/0j9zfsbjor1s6.png)

```wl
Out[]= 0.759528
```

```wl
In[]:= \[ScriptCapitalE] /. \[ScriptCapitalE] -> rval[[nval]]
```

```wl
Out[]= 0.759528
```

#### Level 3

![1njhh5qserhk0](img/1njhh5qserhk0.png)

```wl
Out[]= {0.707106, 0.000486012, 0.00114043, -0.00021351, -0.00021351, 0.000486012, -0.00114043, 0.707106}
```

We need to find the normalization constant for φ(x) , ∞

![1wanj0ojl9w85](img/1wanj0ojl9w85.png)

```wl
Out[]= 4.30598*10^-6
```

![0op12jnnqrv73](img/0op12jnnqrv73.png)

![1kvag7qddsalu](img/1kvag7qddsalu.png)

**Here we can see the bounding state** 

For this first energy level we can see that there is no nodes, since we have two atoms the energy level splits, it is forming a bounding level. This is the origin of the covalent bond which is related with the tunneling effect, since in the forbidden region is non zero, so the electron can be found with probabilities to be in the left region or right region. The blue line is the wave function. 

```wl
In[]:= 
```

![0e4fpgn1ip806](img/0e4fpgn1ip806.png)

```wl
Out[]= 1.
```

What is the probability of finding one electron when $x>a+\frac{b}{2}$?

![1j8w8bz1uk7ae](img/1j8w8bz1uk7ae.png)

```wl
Out[]= 0.387518
```

Let' s find the Δx, then we need the expected value of 〈x〉 and $\left\langle x^2\right\rangle$

![1512pcbbw42x2](img/1512pcbbw42x2.png)

![10niudoykeb3n](img/10niudoykeb3n.png)

```wl
Out[]= 0
```

![1v8wb1l4u9iqa](img/1v8wb1l4u9iqa.png)

```wl
Out[]= 2.73146
```

Then the Δx is :

![13m9pr9kyt2mj](img/13m9pr9kyt2mj.png)

```wl
Out[]= 1.65271
```

Now let' s find the uncertainty in the momentum Δp :

![0hqak8tmkhs66](img/0hqak8tmkhs66.png)

![0s2jt3b00ctzs](img/0s2jt3b00ctzs.png)

```wl
Out[]= 0
```

![13kb9ek884610](img/13kb9ek884610.png)

```wl
Out[]= 3.90633
```

Then Δp is 

![1r70zilb2y59q](img/1r70zilb2y59q.png)

```wl
Out[]= 1.97644
```

Therefore the Heisengerg uncertainty Δx Δp is:

```wl
In[]:= \[CapitalDelta]x \[CapitalDelta]p
```

```wl
Out[]= 3.26649
```

What is the expected value of the kinetic energy k = $\frac{-\hbar ^2}{2m}$

![0ua3s0v5potjs](img/0ua3s0v5potjs.png)

```wl
Out[]= 1.95317
```

What is the expected potential energy for this level?

![1cezxp03uj2w2](img/1cezxp03uj2w2.png)

```wl
Out[]= 0.891571
```

Now, check this out: The total energy of the level

![1705zcawb6ssw](img/1705zcawb6ssw.png)

```wl
Out[]= 2.84474
```

```wl
In[]:= \[ScriptCapitalE] /. \[ScriptCapitalE] -> rval[[nval]]
```

```wl
Out[]= 2.84474
```

#### Level 4

![100tmdl0wefla](img/100tmdl0wefla.png)

```wl
Out[]= {-0.707105, -0.000704019, -0.00116066, 0.000240478, -0.000240478, 0.000704019, -0.00116066, 0.707105}
```

We need to find the normalization constant for φ(x) , ∞

![1igetogowzsae](img/1igetogowzsae.png)

```wl
Out[]= 4.94077*10^-6
```

![1qnh3f3bxi69r](img/1qnh3f3bxi69r.png)

![1ckdmy4sykxat](img/1ckdmy4sykxat.png)

**Here we can see the bounding state** 

For this first energy level we can see that there is no nodes, since we have two atoms the energy level splits, it is forming a bounding level. This is the origin of the covalent bond which is related with the tunneling effect, since in the forbidden region is non zero, so the electron can be found with probabilities to be in the left region or right region. The blue line is the wave function. 

```wl
In[]:= 
```

![18mlahuuiuyg0](img/18mlahuuiuyg0.png)

```wl
Out[]= 1.
```

What is the probability of finding one electron when $x>a+\frac{b}{2}$?

![0dqy6lg2gagzt](img/0dqy6lg2gagzt.png)

```wl
Out[]= 0.345875
```

Let' s find the Δx, then we need the expected value of 〈x〉 and $\left\langle x^2\right\rangle$

![1ispmn5cyrhpc](img/1ispmn5cyrhpc.png)

![0shiztmb8gd3i](img/0shiztmb8gd3i.png)

```wl
Out[]= 0
```

![1vyislb828n2l](img/1vyislb828n2l.png)

```wl
Out[]= 2.87706
```

Then the Δx is :

![17oqhiim3b3zv](img/17oqhiim3b3zv.png)

```wl
Out[]= 1.69619
```

Now let' s find the uncertainty in the momentum Δp :

![0w0quznx6hgeb](img/0w0quznx6hgeb.png)

![1adv3jzn9naw4](img/1adv3jzn9naw4.png)

```wl
Out[]= 0
```

![0wu3jhfacnlqm](img/0wu3jhfacnlqm.png)

```wl
Out[]= 4.58778
```

Then Δp is 

![1ge5ggvgo6ifv](img/1ge5ggvgo6ifv.png)

```wl
Out[]= 2.14191
```

Therefore the Heisengerg uncertainty Δx Δp is:

```wl
In[]:= \[CapitalDelta]x \[CapitalDelta]p
```

```wl
Out[]= 3.63309
```

What is the expected value of the kinetic energy k = $\frac{-\hbar ^2}{2m}$

![0vvbr95gj3vzj](img/0vvbr95gj3vzj.png)

```wl
Out[]= 2.29389
```

What is the expected potential energy for this level?

![04griy39hepah](img/04griy39hepah.png)

```wl
Out[]= 0.670283
```

Now, check this out: The total energy of the level

![0lqu7ac3hsrir](img/0lqu7ac3hsrir.png)

```wl
Out[]= 2.96417
```

```wl
In[]:= \[ScriptCapitalE] /. \[ScriptCapitalE] -> rval[[nval]]
```

```wl
Out[]= 2.96417
```

#### Level 5

![0mgdans32781v](img/0mgdans32781v.png)

```wl
Out[]= {0.705911, -0.0117825, -0.0360802, 0.0157672, 0.0157672, -0.0117825, 0.0360802, 0.705911}
```

We need to find the normalization constant for φ(x) , ∞

![08jhm7qzfpwyu](img/08jhm7qzfpwyu.png)

```wl
Out[]= 0.00528687
```

![1acx3nnir7pga](img/1acx3nnir7pga.png)

![0f93fb1vx41ml](img/0f93fb1vx41ml.png)

**Here we can see the bounding state** 

For this first energy level we can see that there is no nodes, since we have two atoms the energy level splits, it is forming a bounding level. This is the origin of the covalent bond which is related with the tunneling effect, since in the forbidden region is non zero, so the electron can be found with probabilities to be in the left region or right region. The blue line is the wave function. 

```wl
In[]:= 
```

![16bip26iudmvs](img/16bip26iudmvs.png)

```wl
Out[]= 1.
```

What is the probability of finding one electron when $x>a+\frac{b}{2}$?

![1ugxoqr43zp6f](img/1ugxoqr43zp6f.png)

```wl
Out[]= 0.367109
```

Let' s find the Δx, then we need the expected value of 〈x〉 and $\left\langle x^2\right\rangle$

![0hs79k2zoq354](img/0hs79k2zoq354.png)

![1taqkyxy96fja](img/1taqkyxy96fja.png)

```wl
Out[]= 0
```

![0j55ijws4t2ol](img/0j55ijws4t2ol.png)

```wl
Out[]= 3.37312
```

Then the Δx is :

![1c24w0hy53d07](img/1c24w0hy53d07.png)

```wl
Out[]= 1.83661
```

Now let' s find the uncertainty in the momentum Δp :

![1pid7avaeri2m](img/1pid7avaeri2m.png)

![1v5ndohcf8elb](img/1v5ndohcf8elb.png)

```wl
Out[]= 0
```

![15q8ph5amsx36](img/15q8ph5amsx36.png)

```wl
Out[]= 6.17611
```

Then Δp is 

![0ud75w2t8x2bq](img/0ud75w2t8x2bq.png)

```wl
Out[]= 2.48518
```

Therefore the Heisengerg uncertainty Δx Δp is:

```wl
In[]:= \[CapitalDelta]x \[CapitalDelta]p
```

```wl
Out[]= 4.56429
```

What is the expected value of the kinetic energy k = $\frac{-\hbar ^2}{2m}$

![152bpz4z6bn31](img/152bpz4z6bn31.png)

```wl
Out[]= 3.08805
```

What is the expected potential energy for this level?

![1sg833y7dx1xg](img/1sg833y7dx1xg.png)

```wl
Out[]= 2.70104
```

Now, check this out: The total energy of the level

![17hh55s3ht5cp](img/17hh55s3ht5cp.png)

```wl
Out[]= 5.78909
```

```wl
In[]:= \[ScriptCapitalE] /. \[ScriptCapitalE] -> rval[[nval]]
```

```wl
Out[]= 5.78909
```

#### Level 6

![0atxbk1heg2mk](img/0atxbk1heg2mk.png)

```wl
Out[]= {-0.691614, 0.0667736, 0.0750336, -0.10762, 0.10762, -0.0667736, 0.0750336, 0.691614}
```

We need to find the normalization constant for φ(x) , ∞

![16dbz167j7ijr](img/16dbz167j7ijr.png)

```wl
Out[]= 0.0367803
```

![1vzydo1573g8t](img/1vzydo1573g8t.png)

![0ahn1yu3sk8t6](img/0ahn1yu3sk8t6.png)

**Here we can see the bounding state** 

For this first energy level we can see that there is no nodes, since we have two atoms the energy level splits, it is forming a bounding level. This is the origin of the covalent bond which is related with the tunneling effect, since in the forbidden region is non zero, so the electron can be found with probabilities to be in the left region or right region. The blue line is the wave function. 

```wl
In[]:= 
```

![0ik9xa6fn25zx](img/0ik9xa6fn25zx.png)

```wl
Out[]= 1.
```

What is the probability of finding one electron when $x>a+\frac{b}{2}$?

![19od0bijfmxeu](img/19od0bijfmxeu.png)

```wl
Out[]= 0.262557
```

Let' s find the Δx, then we need the expected value of 〈x〉 and $\left\langle x^2\right\rangle$

![0bx70yfsohuc1](img/0bx70yfsohuc1.png)

![0faofkl9sdr4g](img/0faofkl9sdr4g.png)

![1omxqs78l6jfo](img/1omxqs78l6jfo.png)

```wl
Out[]= 0
```

![0q1sp46trklr0](img/0q1sp46trklr0.png)

```wl
Out[]= 4.99442
```

Then the Δx is :

![1r3okk1ia1484](img/1r3okk1ia1484.png)

```wl
Out[]= 2.23482
```

Now let' s find the uncertainty in the momentum Δp :

![184fym8d1fma7](img/184fym8d1fma7.png)

![0mnudk58t3569](img/0mnudk58t3569.png)

```wl
Out[]= 0
```

![0jk8kp1qljc51](img/0jk8kp1qljc51.png)

```wl
Out[]= 7.18132
```

Then Δp is 

![1o6nh12nh92po](img/1o6nh12nh92po.png)

```wl
Out[]= 2.6798
```

Therefore the Heisengerg uncertainty Δx Δp is:

```wl
In[]:= \[CapitalDelta]x \[CapitalDelta]p
```

```wl
Out[]= 5.98887
```

What is the expected value of the kinetic energy k = $\frac{-\hbar ^2}{2m}$

![181vm0xa3vfzn](img/181vm0xa3vfzn.png)

```wl
Out[]= 3.59066
```

What is the expected potential energy for this level?

![1qblg6nr97go9](img/1qblg6nr97go9.png)

```wl
Out[]= 2.60405
```

Now, check this out: The total energy of the level

![1o89jqoobrs56](img/1o89jqoobrs56.png)

```wl
Out[]= 6.19471
```

```wl
In[]:= \[ScriptCapitalE] /. \[ScriptCapitalE] -> rval[[nval]]
```

```wl
Out[]= 6.19471
```

#### Summary of the plots

```wl
In[]:= Show @@ Table[figure[i], {i, 6}]
```

![0cgf4xgi2q8kd](img/0cgf4xgi2q8kd.png)

What is the total energy of the system if it has 6 E?

```wl
In[]:= 
```

What is the total energy of the system if it has 6 E

![1vlmikh3f556y](img/1vlmikh3f556y.png)

```wl
Out[]= 8.68692
```

Kinetic energy of 6 E

```wl
In[]:= 
```

![1wmndjo5dgdiz](img/1wmndjo5dgdiz.png)

```wl
Out[]= 6.21933
```

```wl
In[]:=  
```

Total potential energy of 6 E

![1lkqf8bjr3yjm](img/1lkqf8bjr3yjm.png)

```wl
Out[]= 2.46758
```

#### Cohesion Energy

What is the cohesion energy for this 6 E system,  the cohesion energy is the energy needed to vaporize the system $E_{\text{cohesion} }= \frac{\text{Emol}}{2}-\text{Eatom}$.  What is the formation energy?

![18lgcki9tkdof](img/18lgcki9tkdof.png)

```wl
Out[]= -0.0587472
```

Now what would be the cohesion energy if each atom contributes with 4 E, then the molecule will have 8 E.

![0gmn4xxi34907](img/0gmn4xxi34907.png)

```wl
Out[]= 7.30541
```

![13xgw58doqnd5](img/13xgw58doqnd5.png)

```wl
Out[]= 14.6153
```

![1bylenhb3u6sd](img/1bylenhb3u6sd.png)

```wl
Out[]= 0.00222281
```

Here we notice that  the system evaporates because the cohesion energy is positive meaning that $E_{\text{atom}}$< $E_{\text{mol}}$ 

Now what would be the cohesion energy if each atom contributes with 5 E, then the molecule will have 10 electrons

![1cq9w54hx7tcl](img/1cq9w54hx7tcl.png)

```wl
Out[]= 13.2326
```

![16c2zftfczksm](img/16c2zftfczksm.png)

```wl
Out[]= 26.1935
```

![1waby2wsn7eyt](img/1waby2wsn7eyt.png)

```wl
Out[]= -0.135827
```

```wl
In[]:= 
```

Now we observe that the molecule is more stable. 

```wl
In[]:= 
```

### Three Finite Potential Well: model of a molecule

We want to solve the time independent Schrödinger equation:
$\frac{-\hbar ^2}{2m}$ $\frac{\partial ^2\psi (x)}{\partial x}$ V(x)ψ(x) = E ψ(x)

#### Define the potential

Let's define the finite potential.  use +pw+:

![1afep00t79cgq](img/1afep00t79cgq.png)

![0wkbdax0ymm3b](img/0wkbdax0ymm3b.png)

```wl
In[]:= Plot[v[x] /. {v0 -> vpot, a -> la, b -> lb}, {x, -2 la - lb - 0.3, +2 la + lb + 0.3}, 
   AxesStyle -> {"x", "V(x)"}, 
   PlotStyle -> {Red, Thick}, 
   PlotRange -> {Automatic, {-0.2, vpot + 1.5}}, 
   AxesStyle -> Arrowheads[{0.0, 0.03}], 
   Filling -> Axis 
  ]
```

![0cfttxwgahmg1](img/0cfttxwgahmg1.png)

#### Solving the Schrödinger equation in a . u . 

In atomic units, ℏ -> 1, m -> 1, e -> 1; then the lengths are given in Bohr, 1bohr= 0.53 Å, and the energy in Hartree $E_h$ = 27.21 eV. 
Let's define α = $\sqrt{2(\text{v0} - \mathcal{E})}$, β = $\sqrt{2 \mathcal{E}}$,   (ℰ is +scE+)

![0flc9kg4s1pdt](img/0flc9kg4s1pdt.png)

![000rln0bbersr](img/000rln0bbersr.png)

![0b2y5hzvng7da](img/0b2y5hzvng7da.png)

![1uknqak5xck7b](img/1uknqak5xck7b.png)

![00e6qrp4ve8j7](img/00e6qrp4ve8j7.png)

![03peb6x91039g](img/03peb6x91039g.png)

![17d72gjqh6rm3](img/17d72gjqh6rm3.png)

![04pjld0fio6e9](img/04pjld0fio6e9.png)

![0lhemp1pvohd1](img/0lhemp1pvohd1.png)

Let's define the wavefunction  φ, φ is +j+:

![07eao9rd6t4hb](img/07eao9rd6t4hb.png)

```wl
In[]:= 
```

Now we need to consider the boundary conditions

![1pht867fjlczq](img/1pht867fjlczq.png)

![11ail8revx2il](img/11ail8revx2il.png)

![04n8vtiwms2vk](img/04n8vtiwms2vk.png)

![04il6sqlkhacr](img/04il6sqlkhacr.png)

![1xd9b8o44cfmt](img/1xd9b8o44cfmt.png)

![1unhly4illue7](img/1unhly4illue7.png)

![17t4vfn35ciik](img/17t4vfn35ciik.png)

![0rcg8g17mfsqn](img/0rcg8g17mfsqn.png)

![1lgy5dudqjx26](img/1lgy5dudqjx26.png)

![178peqjb6wtko](img/178peqjb6wtko.png)

![07t7265476wzp](img/07t7265476wzp.png)

![0o3u800j6y14c](img/0o3u800j6y14c.png)

![0xrnd9o6picfg](img/0xrnd9o6picfg.png)

![128s4qssfa97f](img/128s4qssfa97f.png)

```wl
Out[]= {{E^((-a - 2 b) \[Alpha]), -Cos[(-a - 2 b) \[Beta]], -Sin[(-a - 2 b) \[Beta]], 0, 0, 0, 0, 0, 0, 0, 0, 0}, {E^((-a - 2 b) \[Alpha]) \[Alpha], \[Beta] Sin[(-a - 2 b) \[Beta]], -\[Beta] Cos[(-a - 2 b) \[Beta]], 0, 0, 0, 0, 0,0, 0, 0, 0}, {0, Cos[a \[Beta]], -Sin[a \[Beta]], -E^(a \[Alpha]), -E^(-a \[Alpha]), 0, 0, 0, 0, 0, 0, 0}, {0, \[Beta] Sin[a \[Beta]], \[Beta] Cos[a \[Beta]], E^(a \[Alpha]) \[Alpha], -E^(-a \[Alpha]) \[Alpha], 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, E^(b \[Alpha]), E^(-b \[Alpha]), -Cos[b \[Beta]], Sin[b \[Beta]], 0, 0, 0, 0, 0}, {0, 0, 0, -E^(b \[Alpha]) \[Alpha], E^(-b \[Alpha]) \[Alpha], -\[Beta] Sin[b \[Beta]], -\[Beta] Cos[b \[Beta]], 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, Cos[b \[Beta]], Sin[b \[Beta]], -E^(-b \[Alpha]), -E^(b \[Alpha]), 0, 0, 0}, {0, 0, 0, 0, 0, -\[Beta] Sin[b \[Beta]], \[Beta] Cos[b \[Beta]], E^(-b \[Alpha]) \[Alpha], -E^(b \[Alpha]) \[Alpha], 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, E^(-a \[Alpha]), E^(a \[Alpha]), -Cos[a \[Beta]], -Sin[a \[Beta]], 0}, {0, 0, 0, 0, 0, 0, 0, -E^(-a \[Alpha]) \[Alpha], E^(a \[Alpha]) \[Alpha], \[Beta] Sin[a \[Beta]], -\[Beta] Cos[a \[Beta]], 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, Cos[(a + 2 b) \[Beta]],Sin[(a + 2 b) \[Beta]], -E^-((a + 2 b) \[Alpha])}, {0, 0, 0, 0, 0, 0, 0, 0, 0, -\[Beta] Sin[(a + 2 b) \[Beta]], \[Beta] Cos[(a + 2 b) \[Beta]], E^-((a + 2 b) \[Alpha]) \[Alpha]}}
```

```wl
In[]:= mat // MatrixForm
```

|  |  |  |  |  |  |  |  |  |  |  |  |
| - | - | - | - | - | - | - | - | - | - | - | - |
| -Power- | -Times- | -Times- | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| -Times- | -Times- | -Times- | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 0 | -Cos- | -Times- | -Times- | -Times- | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 0 | -Times- | -Times- | -Times- | -Times- | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 0 | 0 | 0 | -Power- | -Power- | -Times- | -Sin- | 0 | 0 | 0 | 0 | 0 |
| 0 | 0 | 0 | -Times- | -Times- | -Times- | -Times- | 0 | 0 | 0 | 0 | 0 |
| 0 | 0 | 0 | 0 | 0 | -Cos- | -Sin- | -Times- | -Times- | 0 | 0 | 0 |
| 0 | 0 | 0 | 0 | 0 | -Times- | -Times- | -Times- | -Times- | 0 | 0 | 0 |
| 0 | 0 | 0 | 0 | 0 | 0 | 0 | -Power- | -Power- | -Times- | -Times- | 0 |
| 0 | 0 | 0 | 0 | 0 | 0 | 0 | -Times- | -Times- | -Times- | -Times- | 0 |
| 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | -Cos- | -Sin- | -Times- |
| 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | -Times- | -Times- | -Times- |

![03mauxto5uxrg](img/03mauxto5uxrg.png)

![0pjqzsusdmj5i](img/0pjqzsusdmj5i.png)

```wl
In[]:= Plot[Det[mat2], {\[ScriptCapitalE], 0, vpot}]
```

![1qb715wdakkud](img/1qb715wdakkud.png)

```wl
In[]:= fig = Plot[Det[mat2], {\[ScriptCapitalE], 0, vpot}, PlotRange -> {-10^-20, 10^-20}]
```

![09r67y8sco37h](img/09r67y8sco37h.png)

```wl
In[]:= fig[[1, 1, 1, 1, 1, 3, 1]][[2 ;; -1]]
```

```wl
Out[]= {Line[{{5.73363, 1.*10^-20}, {5.73363, -1.*10^-20}}], Line[{{5.98363, -1.*10^-20}, {5.98363, 1.*10^-20}}], Line[{{2.90276, 1.*10^-20}, {2.90276, -1.*10^-20}}], Line[{{0.763652, 1.*10^-20}, {0.763652, -1.*10^-20}}], Line[{{6.34613, 1.*10^-20}, {6.34613, -1.*10^-20}}], Line[{{2.82144, -1.*10^-20}, {2.82144, 1.*10^-20}}], Line[{{2.99045, -1.*10^-20}, {2.99045, 1.*10^-20}}], Line[{{0.735003, 1.*10^-20}, {0.735003, -1.*10^-20}}], Line[{{0.749238, -1.*10^-20}, {0.749238, 1.*10^-20}}]}
```

```wl
In[]:= approxroot = fig[[1, 1, 1, 1, 1, 3, 1]][[2 ;; -1]] // Map[#[[1, 1, 1]] &, #] & // Sort
```

```wl
Out[]= {0.735003, 0.749238, 0.763652, 2.82144, 2.90276, 2.99045, 5.73363, 5.98363, 6.34613}
```

```wl
In[]:= rval = approxroot // Map[FindRoot[Det[mat2], {\[ScriptCapitalE], #}] &, #] & // Map[\[ScriptCapitalE] /. # &, #] &
```

```wl
Out[]= {0.734976, 0.749237, 0.76373, 2.82139, 2.90277, 2.99047, 5.73363, 5.98363, 6.34616}
```

```wl
In[]:= fig[[1, 1, 1]]
```

![0gprxnchp7hje](img/0gprxnchp7hje.png)

#### Level 1

![04yvqzpubnp0s](img/04yvqzpubnp0s.png)

![18zeg3up5kans](img/18zeg3up5kans.png)

![1ccdhirc4q0y2](img/1ccdhirc4q0y2.png)

![0rfuuk1dzpu1z](img/0rfuuk1dzpu1z.png)

![0amlm9rhegdsf](img/0amlm9rhegdsf.png)

```wl
Out[]= {0.707107, -2.35466*10^-6, 1.22653*10^-6, 1.04049*10^-9, 0.0000387152, 3.78837*10^-6, 0, 0.0000387152, 1.04049*10^-9, -2.35466*10^-6, -1.22652*10^-6, 0.707106}
```

We need to find the normalization constant for φ(x) , ∞

![0s5py3943a861](img/0s5py3943a861.png)

```wl
Out[]= 3.72987*10^-11
```

![0fr96c5mjrj6b](img/0fr96c5mjrj6b.png)

```wl
Out[]= 163739.
```

![159q3204u0t0r](img/159q3204u0t0r.png)

![1ib0k9976dfwx](img/1ib0k9976dfwx.png)

**Here we can see the bounding state** 

For this first energy level we can see that there is no nodes, since we have two atoms the energy level splits, it is forming a bounding level. This is the origin of the covalent bond which is related with the tunneling effect, since in the forbidden region is non zero, so the electron can be found with probabilities to be in the left region or right region. The blue line is the wave function. 

```wl
In[]:= 
```

![0ma9fkicd2q6b](img/0ma9fkicd2q6b.png)

```wl
Out[]= 1.
```

What is the probability of finding one electron when $x>a+\frac{b}{2}$?

![1vpkirywvafon](img/1vpkirywvafon.png)

```wl
Out[]= 0.489026
```

Let' s find the Δx, then we need the expected value of 〈x〉 and $\left\langle x^2\right\rangle$

![01n4srwdmzn3c](img/01n4srwdmzn3c.png)

```wl
Out[]= -2.60588*10^-6
```

![1cco1712whwo8](img/1cco1712whwo8.png)

```wl
Out[]= 4.64783
```

Then the Δx is :

![13fjmatmgqg1c](img/13fjmatmgqg1c.png)

```wl
Out[]= 2.15588
```

Now let' s find the uncertainty in the momentum Δp :

![0u9k8dc2edchd](img/0u9k8dc2edchd.png)

```wl
Out[]= 0. - 6.74607*10^-9 I
```

![0vef80kscf1t0](img/0vef80kscf1t0.png)

```wl
Out[]= 1.07139
```

Then Δp is 

![0owvcnx915lzl](img/0owvcnx915lzl.png)

```wl
Out[]= 1.03508 + 0. I
```

Therefore the Heisengerg uncertainty Δx Δp is:

```wl
In[]:= \[CapitalDelta]x \[CapitalDelta]p
```

```wl
Out[]= 2.23151 + 0. I
```

What is the expected value of the kinetic energy k = $\frac{-\hbar ^2}{2m}$

![07ipb1h349fva](img/07ipb1h349fva.png)

```wl
Out[]= 0.535694
```

What is the expected potential energy for this level?

![0v4lgsnc8assd](img/0v4lgsnc8assd.png)

```wl
Out[]= 0.199283
```

Now, check this out: The total energy of the level

![06o73zl5yndj5](img/06o73zl5yndj5.png)

```wl
Out[]= 0.734976
```

```wl
In[]:= \[ScriptCapitalE] /. \[ScriptCapitalE] -> rval[[nval]]
```

```wl
Out[]= 0.734976
```

#### Level 2

![0acrervd5zcsl](img/0acrervd5zcsl.png)

```wl
Out[]= {-0.707107, 2.3066*10^-6, -1.35301*10^-6, -1.02945*10^-9, -6.98887*10^-7, 0, 5.75311*10^-8, 6.9892*10^-7, 1.02945*10^-9, -2.3066*10^-6, -1.35301*10^-6, 0.707107}
```

We need to find the normalization constant for φ(x) , ∞

![13yx78d7gxzzi](img/13yx78d7gxzzi.png)

```wl
Out[]= 1.85273*10^-11
```

![0us88buhkwy8n](img/0us88buhkwy8n.png)

```wl
Out[]= 232324.
```

![1qqz73p2bd3sz](img/1qqz73p2bd3sz.png)

![0u0sq9shoyn4p](img/0u0sq9shoyn4p.png)

**Here we can see the bounding state** 

For this first energy level we can see that there is no nodes, since we have two atoms the energy level splits, it is forming a bounding level. This is the origin of the covalent bond which is related with the tunneling effect, since in the forbidden region is non zero, so the electron can be found with probabilities to be in the left region or right region. The blue line is the wave function. 

```wl
In[]:= 
```

![19snkjry4ea5i](img/19snkjry4ea5i.png)

```wl
Out[]= 1.
```

What is the probability of finding one electron when $x>a+\frac{b}{2}$?

![1dkw740yak5w8](img/1dkw740yak5w8.png)

```wl
Out[]= 0.00013201
```

Let' s find the Δx, then we need the expected value of 〈x〉 and $\left\langle x^2\right\rangle$

![07r61podp7wdv](img/07r61podp7wdv.png)

```wl
Out[]= 3.07691*10^-7
```

![149a8ri38olo3](img/149a8ri38olo3.png)

```wl
Out[]= 9.22532
```

Then the Δx is :

![1auso8ewyunef](img/1auso8ewyunef.png)

```wl
Out[]= 3.03732
```

Now let' s find the uncertainty in the momentum Δp :

![0o8xc483l68rx](img/0o8xc483l68rx.png)

![0uidco8imwv23](img/0uidco8imwv23.png)

```wl
Out[]= 0
```

![14sgdkmt4ggrm](img/14sgdkmt4ggrm.png)

```wl
Out[]= 1.15522
```

Then Δp is 

![0sp64d66nh2dd](img/0sp64d66nh2dd.png)

```wl
Out[]= 1.07481
```

Therefore the Heisengerg uncertainty Δx Δp is:

```wl
In[]:= \[CapitalDelta]x \[CapitalDelta]p
```

```wl
Out[]= 3.26455
```

What is the expected value of the kinetic energy k = $\frac{-\hbar ^2}{2m}$

![1o54oq3pt5rac](img/1o54oq3pt5rac.png)

```wl
Out[]= 0.577611
```

What is the expected potential energy for this level?

![1xkjml3zkfou9](img/1xkjml3zkfou9.png)

```wl
Out[]= 0.171626
```

Now, check this out: The total energy of the level

![0mhs651syyjrc](img/0mhs651syyjrc.png)

```wl
Out[]= 0.749237
```

```wl
In[]:= \[ScriptCapitalE] /. \[ScriptCapitalE] -> rval[[nval]]
```

```wl
Out[]= 0.749237
```

#### Level 3

![07zwy5ngf932e](img/07zwy5ngf932e.png)

```wl
Out[]= {0.707107, -2.25209*10^-6, 1.47904*10^-6, 1.01746*10^-9, -0.0000375897, -3.77483*10^-6, 0, -0.0000375897, 1.01746*10^-9, -2.25209*10^-6, -1.47904*10^-6, 0.707106}
```

We need to find the normalization constant for φ(x) , ∞

![0ewfq73vv9n1i](img/0ewfq73vv9n1i.png)

```wl
Out[]= 3.68009*10^-11
```

![0dztrjul5wylp](img/0dztrjul5wylp.png)

```wl
Out[]= 164843.
```

![1wrm0gjpszro9](img/1wrm0gjpszro9.png)

![0762xeubg30vq](img/0762xeubg30vq.png)

**Here we can see the bounding state** 

For this first energy level we can see that there is no nodes, since we have two atoms the energy level splits, it is forming a bounding level. This is the origin of the covalent bond which is related with the tunneling effect, since in the forbidden region is non zero, so the electron can be found with probabilities to be in the left region or right region. The blue line is the wave function. 

```wl
In[]:= 
```

![1nybsm4tsv6k2](img/1nybsm4tsv6k2.png)

```wl
Out[]= 1.
```

What is the probability of finding one electron when $x>a+\frac{b}{2}$?

![16sjw2whmsoze](img/16sjw2whmsoze.png)

```wl
Out[]= 0.484451
```

Let' s find the Δx, then we need the expected value of 〈x〉 and $\left\langle x^2\right\rangle$

![05jxdgt565f8u](img/05jxdgt565f8u.png)

```wl
Out[]= -1.914*10^-6
```

![16xkbd7kjjpxv](img/16xkbd7kjjpxv.png)

```wl
Out[]= 4.81129
```

Then the Δx is :

![10xks8rlh1eat](img/10xks8rlh1eat.png)

```wl
Out[]= 2.19346
```

Now let' s find the uncertainty in the momentum Δp :

![1xkui8n61dbms](img/1xkui8n61dbms.png)

```wl
Out[]= 0. + 3.89293*10^-9 I
```

![1sd1zt1ww8a71](img/1sd1zt1ww8a71.png)

```wl
Out[]= 1.24175
```

Then Δp is 

![117r6agry8ih2](img/117r6agry8ih2.png)

```wl
Out[]= 1.11434 + 0. I
```

Therefore the Heisengerg uncertainty Δx Δp is:

```wl
In[]:= \[CapitalDelta]x \[CapitalDelta]p
```

```wl
Out[]= 2.44426 + 0. I
```

What is the expected value of the kinetic energy k = $\frac{-\hbar ^2}{2m}$

![0cjo3ikzmqnph](img/0cjo3ikzmqnph.png)

```wl
Out[]= 0.620876
```

What is the expected potential energy for this level?

![10xrzwib3ge4i](img/10xrzwib3ge4i.png)

```wl
Out[]= 0.142854
```

Now, check this out: The total energy of the level

![1cy0z55o8q1fs](img/1cy0z55o8q1fs.png)

```wl
Out[]= 0.76373
```

```wl
In[]:= \[ScriptCapitalE] /. \[ScriptCapitalE] -> rval[[nval]]
```

```wl
Out[]= 0.76373
```

#### Level 4

![1mi8ibyo0en6b](img/1mi8ibyo0en6b.png)

```wl
Out[]= {-0.707107, 0.0000148909, 0.0000145698, 6.09589*10^-8, 0.000294205, 0, -0.0000294891, -0.000294205, -6.09589*10^-8, -0.0000148909, 0.0000145698, 0.707107}
```

We need to find the normalization constant for φ(x) , ∞

![0ozlcnq072on7](img/0ozlcnq072on7.png)

```wl
Out[]= 2.45375*10^-9
```

![12vl3lkq2dca5](img/12vl3lkq2dca5.png)

```wl
Out[]= 20187.6
```

![0vgqar98rag0v](img/0vgqar98rag0v.png)

![0uf1g62qmqspa](img/0uf1g62qmqspa.png)

**Here we can see the bounding state** 

For this first energy level we can see that there is no nodes, since we have two atoms the energy level splits, it is forming a bounding level. This is the origin of the covalent bond which is related with the tunneling effect, since in the forbidden region is non zero, so the electron can be found with probabilities to be in the left region or right region. The blue line is the wave function. 

```wl
In[]:= 
```

![0wki72dct73km](img/0wki72dct73km.png)

```wl
Out[]= 1.
```

What is the probability of finding one electron when $x>a+\frac{b}{2}$?

![1bp3osyztdypn](img/1bp3osyztdypn.png)

```wl
Out[]= 0.42894
```

Let' s find the Δx, then we need the expected value of 〈x〉 and $\left\langle x^2\right\rangle$

![1risa3vko0k50](img/1risa3vko0k50.png)

```wl
Out[]= 2.58958*10^-8
```

![1qzexv45wi8my](img/1qzexv45wi8my.png)

```wl
Out[]= 4.9508
```

Then the Δx is :

![1c3jfcxkwjccm](img/1c3jfcxkwjccm.png)

```wl
Out[]= 2.22504
```

Now let' s find the uncertainty in the momentum Δp :

![1twh6xzw14n4o](img/1twh6xzw14n4o.png)

```wl
Out[]= 0. + 2.68603*10^-10 I
```

![1dk7sphc69r4j](img/1dk7sphc69r4j.png)

```wl
Out[]= 3.78088
```

Then Δp is 

![1sjpzsgw0qpvk](img/1sjpzsgw0qpvk.png)

```wl
Out[]= 1.94445 + 0. I
```

Therefore the Heisengerg uncertainty Δx Δp is:

```wl
In[]:= \[CapitalDelta]x \[CapitalDelta]p
```

```wl
Out[]= 4.32647 + 0. I
```

What is the expected value of the kinetic energy k = $\frac{-\hbar ^2}{2m}$

![1gohuio9x1wl5](img/1gohuio9x1wl5.png)

```wl
Out[]= 1.89044
```

What is the expected potential energy for this level?

![1dznjp1hbb11c](img/1dznjp1hbb11c.png)

```wl
Out[]= 0.93095
```

Now, check this out: The total energy of the level

![0wk4astqskkxo](img/0wk4astqskkxo.png)

```wl
Out[]= 2.82139
```

```wl
In[]:= \[ScriptCapitalE] /. \[ScriptCapitalE] -> rval[[nval]]
```

```wl
Out[]= 2.82139
```

#### Level 5

![1n02rozg5qnx3](img/1n02rozg5qnx3.png)

```wl
Out[]= {0.707107, -0.0000187817, -0.0000135742, -7.24756*10^-8, -1.6508*10^-6, 1.57629*10^-6, 0, -1.6508*10^-6, -7.24756*10^-8, -0.0000187817, 0.0000135742, 0.707107}
```

We need to find the normalization constant for φ(x) , ∞

![03cb7cakuyjop](img/03cb7cakuyjop.png)

```wl
Out[]= 1.4765*10^-9
```

![1k0dhec2rgb3v](img/1k0dhec2rgb3v.png)

```wl
Out[]= 26024.6
```

![1oltf7r19mx0d](img/1oltf7r19mx0d.png)

![1lkfoi8npuuqz](img/1lkfoi8npuuqz.png)

**Here we can see the bounding state** 

For this first energy level we can see that there is no nodes, since we have two atoms the energy level splits, it is forming a bounding level. This is the origin of the covalent bond which is related with the tunneling effect, since in the forbidden region is non zero, so the electron can be found with probabilities to be in the left region or right region. The blue line is the wave function. 

```wl
In[]:= 
```

![1dt2xil97u6nv](img/1dt2xil97u6nv.png)

```wl
Out[]= 1.
```

What is the probability of finding one electron when $x>a+\frac{b}{2}$?

![07wt4ruuvax1s](img/07wt4ruuvax1s.png)

```wl
Out[]= 0.00133559
```

Let' s find the Δx, then we need the expected value of 〈x〉 and $\left\langle x^2\right\rangle$

![1xcmsi3gxn9pi](img/1xcmsi3gxn9pi.png)

```wl
Out[]= -2.69909*10^-9
```

![14qtpg457jpfh](img/14qtpg457jpfh.png)

```wl
Out[]= 9.53588
```

Then the Δx is :

![0a9ngeyyy78eg](img/0a9ngeyyy78eg.png)

```wl
Out[]= 3.08802
```

Now let' s find the uncertainty in the momentum Δp :

![0jtvbdlkdipwm](img/0jtvbdlkdipwm.png)

![07cvrlgrd51z6](img/07cvrlgrd51z6.png)

```wl
Out[]= 0
```

![1i7wr9hhuc1az](img/1i7wr9hhuc1az.png)

```wl
Out[]= 4.23045
```

Then Δp is 

![1jyy6fqn6xk4q](img/1jyy6fqn6xk4q.png)

```wl
Out[]= 2.05681
```

Therefore the Heisengerg uncertainty Δx Δp is:

```wl
In[]:= \[CapitalDelta]x \[CapitalDelta]p
```

```wl
Out[]= 6.35147
```

What is the expected value of the kinetic energy k = $\frac{-\hbar ^2}{2m}$

![1ccqbe9kmz2h6](img/1ccqbe9kmz2h6.png)

```wl
Out[]= 2.11523
```

What is the expected potential energy for this level?

![07l0hevdvzf6m](img/07l0hevdvzf6m.png)

```wl
Out[]= 0.787543
```

Now, check this out: The total energy of the level

![0ean2n18ialky](img/0ean2n18ialky.png)

```wl
Out[]= 2.90277
```

```wl
In[]:= \[ScriptCapitalE] /. \[ScriptCapitalE] -> rval[[nval]]
```

```wl
Out[]= 2.90277
```

#### Level 6

![1lre4ahv4imsg](img/1lre4ahv4imsg.png)

```wl
Out[]= {-0.707107, 0.0000232449, 0.0000117421, 8.71578*10^-8, -0.000349959, 0, 0.0000366659, 0.000349959, -8.71578*10^-8, -0.0000232449, 0.0000117421, 0.707107}
```

We need to find the normalization constant for φ(x) , ∞

![1vmlliapurybu](img/1vmlliapurybu.png)

```wl
Out[]= 3.58269*10^-9
```

![13hexb7wpa36q](img/13hexb7wpa36q.png)

```wl
Out[]= 16706.9
```

![0wc369cvcr2vs](img/0wc369cvcr2vs.png)

![1d9vu5q5lc8ew](img/1d9vu5q5lc8ew.png)

**Here we can see the bounding state** 

For this first energy level we can see that there is no nodes, since we have two atoms the energy level splits, it is forming a bounding level. This is the origin of the covalent bond which is related with the tunneling effect, since in the forbidden region is non zero, so the electron can be found with probabilities to be in the left region or right region. The blue line is the wave function. 

```wl
In[]:= 
```

![0zuw82k389jrz](img/0zuw82k389jrz.png)

```wl
Out[]= 1.
```

What is the probability of finding one electron when $x>a+\frac{b}{2}$?

![0kzwj8e0ir64u](img/0kzwj8e0ir64u.png)

```wl
Out[]= 0.450739
```

Let' s find the Δx, then we need the expected value of 〈x〉 and $\left\langle x^2\right\rangle$

![1s3g63ea0ci0k](img/1s3g63ea0ci0k.png)

```wl
Out[]= -5.26773*10^-9
```

![1mnma5qlicshw](img/1mnma5qlicshw.png)

```wl
Out[]= 5.18218
```

Then the Δx is :

![19w27c0kyqx2n](img/19w27c0kyqx2n.png)

```wl
Out[]= 2.27644
```

Now let' s find the uncertainty in the momentum Δp :

![19phaecvowbbu](img/19phaecvowbbu.png)

```wl
Out[]= 0. - 1.4165*10^-10 I
```

![0ref2jqh4jcsl](img/0ref2jqh4jcsl.png)

```wl
Out[]= 4.74773
```

Then Δp is 

![0j8dwcp8tqqm2](img/0j8dwcp8tqqm2.png)

```wl
Out[]= 2.17893 + 0. I
```

Therefore the Heisengerg uncertainty Δx Δp is:

```wl
In[]:= \[CapitalDelta]x \[CapitalDelta]p
```

```wl
Out[]= 4.9602 + 0. I
```

What is the expected value of the kinetic energy k = $\frac{-\hbar ^2}{2m}$

![0asdxrksackfu](img/0asdxrksackfu.png)

```wl
Out[]= 2.37386
```

What is the expected potential energy for this level?

![0f1mfxmnhe9o4](img/0f1mfxmnhe9o4.png)

```wl
Out[]= 0.616609
```

Now, check this out: The total energy of the level

![0igqrk5l1pvsm](img/0igqrk5l1pvsm.png)

```wl
Out[]= 2.99047
```

```wl
In[]:= \[ScriptCapitalE] /. \[ScriptCapitalE] -> rval[[nval]]
```

```wl
Out[]= 2.99047
```

#### Level 7

![11485a9qnmka3](img/11485a9qnmka3.png)

```wl
Out[]= {0.706805, 0.00430286, -0.00312713, 0.000306164, 0.0193428, -0.00686909, 0, 0.0193428, 0.000306164, 0.00430286, 0.00312713, 0.706805}
```

We need to find the normalization constant for φ(x) , ∞

![0tnm5n821l4ie](img/0tnm5n821l4ie.png)

```wl
Out[]= 0.000187818
```

![0vz8g31njb006](img/0vz8g31njb006.png)

```wl
Out[]= 72.9679
```

![0grnprpmudlb2](img/0grnprpmudlb2.png)

![0nqlftq82wrc8](img/0nqlftq82wrc8.png)

**Here we can see the bounding state** 

For this first energy level we can see that there is no nodes, since we have two atoms the energy level splits, it is forming a bounding level. This is the origin of the covalent bond which is related with the tunneling effect, since in the forbidden region is non zero, so the electron can be found with probabilities to be in the left region or right region. The blue line is the wave function. 

```wl
In[]:= 
```

![13lhrpfwnj61u](img/13lhrpfwnj61u.png)

```wl
Out[]= 1.
```

What is the probability of finding one electron when $x>a+\frac{b}{2}$?

![1vzvmipqjgha8](img/1vzvmipqjgha8.png)

```wl
Out[]= 0.268665
```

Let' s find the Δx, then we need the expected value of 〈x〉 and $\left\langle x^2\right\rangle$

![0e1zdy8xlvabu](img/0e1zdy8xlvabu.png)

![1iwzjtwm4vui9](img/1iwzjtwm4vui9.png)

```wl
Out[]= 0
```

![17c92kaum82ju](img/17c92kaum82ju.png)

```wl
Out[]= 5.94924
```

Then the Δx is :

![17fjka5v8c5be](img/17fjka5v8c5be.png)

```wl
Out[]= 2.43911
```

Now let' s find the uncertainty in the momentum Δp :

![09lzlk08mvktc](img/09lzlk08mvktc.png)

![0jvi1jycq3ohj](img/0jvi1jycq3ohj.png)

```wl
Out[]= 0
```

![030sl1ft55bhp](img/030sl1ft55bhp.png)

```wl
Out[]= 6.14248
```

Then Δp is 

![0jxqbjnqnkovx](img/0jxqbjnqnkovx.png)

```wl
Out[]= 2.4784
```

Therefore the Heisengerg uncertainty Δx Δp is:

```wl
In[]:= \[CapitalDelta]x \[CapitalDelta]p
```

```wl
Out[]= 6.04509
```

What is the expected value of the kinetic energy k = $\frac{-\hbar ^2}{2m}$

![133nkbgb8h7gp](img/133nkbgb8h7gp.png)

```wl
Out[]= 3.07124
```

What is the expected potential energy for this level?

![1ck3oiuuwf79d](img/1ck3oiuuwf79d.png)

```wl
Out[]= 2.66239
```

Now, check this out: The total energy of the level

![00ufldwjb7g4s](img/00ufldwjb7g4s.png)

```wl
Out[]= 5.73363
```

```wl
In[]:= \[ScriptCapitalE] /. \[ScriptCapitalE] -> rval[[nval]]
```

```wl
Out[]= 5.73363
```

#### Level 8

![11if7358aizc5](img/11if7358aizc5.png)

```wl
Out[]= {-0.706904, -0.00698728, 0.0105403, -0.0017447, 0.0109703, 0, -0.00271858, -0.0109703, 0.0017447, 0.00698728, 0.0105403, 0.706904}
```

We need to find the normalization constant for φ(x) , ∞

![0jfu6fqkpon45](img/0jfu6fqkpon45.png)

```wl
Out[]= 0.000587414
```

![15j4ovhbg717x](img/15j4ovhbg717x.png)

```wl
Out[]= 41.2599
```

![1sj1o58tb6n69](img/1sj1o58tb6n69.png)

![0d1hxrr9g0rf5](img/0d1hxrr9g0rf5.png)

**Here we can see the bounding state** 

For this first energy level we can see that there is no nodes, since we have two atoms the energy level splits, it is forming a bounding level. This is the origin of the covalent bond which is related with the tunneling effect, since in the forbidden region is non zero, so the electron can be found with probabilities to be in the left region or right region. The blue line is the wave function. 

```wl
In[]:= 
```

![03n8xkbbp6mjm](img/03n8xkbbp6mjm.png)

```wl
Out[]= 1.
```

What is the probability of finding one electron when $x>a+\frac{b}{2}$?

![13ii68pc8sfqf](img/13ii68pc8sfqf.png)

```wl
Out[]= 0.0115022
```

Let' s find the Δx, then we need the expected value of 〈x〉 and $\left\langle x^2\right\rangle$

![0lgbtznm6db02](img/0lgbtznm6db02.png)

![1s5t9b3dckuiu](img/1s5t9b3dckuiu.png)

```wl
Out[]= 0
```

![1az3cilpb31hj](img/1az3cilpb31hj.png)

```wl
Out[]= 11.048
```

Then the Δx is :

![0jk7s3sm1lmal](img/0jk7s3sm1lmal.png)

```wl
Out[]= 3.32386
```

Now let' s find the uncertainty in the momentum Δp :

![03i292e5hbrpw](img/03i292e5hbrpw.png)

![11gopas8hiqod](img/11gopas8hiqod.png)

```wl
Out[]= 0
```

![1tgm7ypapl7ss](img/1tgm7ypapl7ss.png)

```wl
Out[]= 6.80117
```

Then Δp is 

![130rtyj9ymj71](img/130rtyj9ymj71.png)

```wl
Out[]= 2.60791
```

Therefore the Heisengerg uncertainty Δx Δp is:

```wl
In[]:= \[CapitalDelta]x \[CapitalDelta]p
```

```wl
Out[]= 8.66831
```

What is the expected value of the kinetic energy k = $\frac{-\hbar ^2}{2m}$

![00ooh3onxa0aa](img/00ooh3onxa0aa.png)

```wl
Out[]= 3.40059
```

What is the expected potential energy for this level?

![0g4lsmdplzvkm](img/0g4lsmdplzvkm.png)

```wl
Out[]= 2.58307
```

Now, check this out: The total energy of the level

![18a3zjra365u2](img/18a3zjra365u2.png)

```wl
Out[]= 5.98366
```

```wl
In[]:= \[ScriptCapitalE] /. \[ScriptCapitalE] -> rval[[nval]]
```

```wl
Out[]= 5.98363
```

#### Level 9

![1x78aq9incjz4](img/1x78aq9incjz4.png)

```wl
Out[]= {0.610895, 0.00276313, -0.067175, 0.0537556, -0.33678, 0.109335, 0, -0.33678, 0.0537556, 0.00276313, 0.067175, 0.610895}
```

We need to find the normalization constant for φ(x) , ∞

![1vg8s3uyvaxsy](img/1vg8s3uyvaxsy.png)

```wl
Out[]= 0.0359037
```

![14gc63qgj4kq7](img/14gc63qgj4kq7.png)

```wl
Out[]= 5.27753
```

![13zvamkuwbi5i](img/13zvamkuwbi5i.png)

![0iaw88mihf3vd](img/0iaw88mihf3vd.png)

**Here we can see the bounding state** 

For this first energy level we can see that there is no nodes, since we have two atoms the energy level splits, it is forming a bounding level. This is the origin of the covalent bond which is related with the tunneling effect, since in the forbidden region is non zero, so the electron can be found with probabilities to be in the left region or right region. The blue line is the wave function. 

```wl
In[]:= 
```

![1cv5sjtj6yn6l](img/1cv5sjtj6yn6l.png)

```wl
Out[]= 1.
```

What is the probability of finding one electron when $x>a+\frac{b}{2}$?

![1bt54jwx93vjc](img/1bt54jwx93vjc.png)

```wl
Out[]= 0.36781
```

Let' s find the Δx, then we need the expected value of 〈x〉 and $\left\langle x^2\right\rangle$

![1u19vno3yar5s](img/1u19vno3yar5s.png)

![0hqk6zh2stact](img/0hqk6zh2stact.png)

```wl
Out[]= 0
```

![036gobw0jnwo7](img/036gobw0jnwo7.png)

```wl
Out[]= 8.40386
```

Then the Δx is :

![058zwogq4z8w6](img/058zwogq4z8w6.png)

```wl
Out[]= 2.89894
```

Now let' s find the uncertainty in the momentum Δp :

![1u2o9tmjbbid1](img/1u2o9tmjbbid1.png)

![060p6n70gpj9p](img/060p6n70gpj9p.png)

```wl
Out[]= 0
```

![0jbhfpwwui7mi](img/0jbhfpwwui7mi.png)

```wl
Out[]= 8.04229
```

Then Δp is 

![1qkalppttts1b](img/1qkalppttts1b.png)

```wl
Out[]= 2.83589
```

Therefore the Heisengerg uncertainty Δx Δp is:

```wl
In[]:= \[CapitalDelta]x \[CapitalDelta]p
```

```wl
Out[]= 8.22109
```

What is the expected value of the kinetic energy k = $\frac{-\hbar ^2}{2m}$

![042ccge4e5056](img/042ccge4e5056.png)

```wl
Out[]= 4.02115
```

What is the expected potential energy for this level?

![097ibi8xm6pou](img/097ibi8xm6pou.png)

```wl
Out[]= 2.32501
```

Now, check this out: The total energy of the level

![1a54iu84gveuk](img/1a54iu84gveuk.png)

```wl
Out[]= 6.34616
```

```wl
In[]:= \[ScriptCapitalE] /. \[ScriptCapitalE] -> rval[[nval]]
```

```wl
Out[]= 6.34616
```

#### Summary of the plots

```wl
In[]:= Show @@ Table[figure[i], {i, 9}]
```

![19k5la0cfk93r](img/19k5la0cfk93r.png)

What is the total energy of the system if it has 6 E?

```wl
In[]:= 
```

What is the total energy of the system if it has 6 E

![1vlmikh3f556y](img/1vlmikh3f556y.png)

```wl
Out[]= 8.68692
```

Kinetic energy of 6 E

```wl
In[]:= 
```

![1wmndjo5dgdiz](img/1wmndjo5dgdiz.png)

```wl
Out[]= 6.21933
```

```wl
In[]:=  
```

Total potential energy of 6 E

![1lkqf8bjr3yjm](img/1lkqf8bjr3yjm.png)

```wl
Out[]= 2.46758
```

#### Cohesion Energy

What is the cohesion energy for this 6 E system,  the cohesion energy is the energy needed to vaporize the system $E_{\text{cohesion} }= \frac{\text{Emol}}{2}-\text{Eatom}$.  What is the formation energy?

![18lgcki9tkdof](img/18lgcki9tkdof.png)

```wl
Out[]= -0.0587472
```

Now what would be the cohesion energy if each atom contributes with 4 E, then the molecule will have 8 E.

![0gmn4xxi34907](img/0gmn4xxi34907.png)

```wl
Out[]= 7.30541
```

![13xgw58doqnd5](img/13xgw58doqnd5.png)

```wl
Out[]= 14.6153
```

![1bylenhb3u6sd](img/1bylenhb3u6sd.png)

```wl
Out[]= 0.00222281
```

Here we notice that  the system evaporates because the cohesion energy is positive meaning that $E_{\text{atom}}$< $E_{\text{mol}}$ 

Now what would be the cohesion energy if each atom contributes with 5 E, then the molecule will have 10 electrons

![1cq9w54hx7tcl](img/1cq9w54hx7tcl.png)

```wl
Out[]= 13.2326
```

![16c2zftfczksm](img/16c2zftfczksm.png)

```wl
Out[]= 26.1935
```

![1waby2wsn7eyt](img/1waby2wsn7eyt.png)

```wl
Out[]= -0.135827
```

```wl
In[]:= 
```

Now we observe that the molecule is more stable. 

```wl
In[]:= 
```