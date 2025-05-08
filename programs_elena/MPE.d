'H42b','H42a','H42c'     ! (f6) – names of history files (average values of energy, order parameters, temperature etc.)

17.001 -2.0 8   ! To, dT, nT (initial temperature, step and number of temp. steps)

'r42','r42','r42','r42','r42','r42','r42','r42','r42','r42','r42','r42','r42','r42','r42','r42','r42','r42','r42','r42','r42','r42','r42','r42','r42','r42','r42','r42','r42','r42','r42','r42','r42','r42','r42','r42','r42','r42','r42','r42','r42','r42','r42','r42','r42','r42','r42','r42','r42','r42' ! (configuration files for each temperature, each file has format x(i),y(i),z(i),alpha(i), beta(i), gamma(i), energy(i) for each atomic site)

'prov43'         ! s4     - temporary configuration for any case 
 
'r'            ! (f1) – technical data not important

1              ! (Nj) - Anzahl von Zeilen in J (how many different l-s (Q_lm) I have; For example I have at each site a complex multipole consisting of Q_30, Q_3-1 and Q_20, in that case Nj=2)

'4'            ! ?? (f2) - L, which ls I have (in the above example I have 3,2)

1             ! (Njm)- Anzahl von Paaren in JM (how many different lm combinations (Q_lm) I have; in the example above Njm=3)

'4 0 1'  ! (f3) - LM pairs explicitly

'AmmBen200'        ! (f4) – input file with coordinates and spatial orientations of moments (format x(i),y(i),z(i),alpha(i), beta(i), gamma(i), energy(i) for each atomic site)

0           ! Hin – external magnetic field in plane

0           ! angle between Ox and in-plane field

0           !Hz – out-of=plane field

257            ! (Ns) - Anzahl von Momenten in NM

500           ! (nmcs)-Numb.of MC steps

500            ! (nmcr)-numb.of MC steps for writing

8             ! (init)-nombre aleatoire initial (initial random number)

0          ! (iocSM ig) - iocSM ig=0, aleat| iconfig=1, .init.donnee' (0 – random initial configuration, 1 – certain configuration)

