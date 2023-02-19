sets
         i locations /1*13/
         t Days      /1*5/;
alias
         (i,j), (t,k,l);
scalar
         C  Transportation cost  /5400/
         Ci Inventory cost per day /10/
         Cs Shortage cost per each shortage per distance /5/    ;
parameter
         d1(i) Distances between regions and warehouse
/
1        0
2        155
3        105
4        105
5        155
6        115
7        65
8        65
9        115
10       75
11       25
12       25
13       75
/;
table  dprime(i,t)   the real demand
       1     2     3     4     5
1      0     0     0     0     0
2      108   91    113   124   101
3      21    20    18    23    22
4      55    61    58    54    66
5      62    67    67    58    62
6      60    58    67    61    52
7      61    49    68    61    56
8      37    40    49    48    32
9      42    34    30    37    35
10     57    62    60    60    60
11     84    78    91    85    85
12     73    76    84    85    82
13     55    72    65    72    52

;
table  d(i,t)       the predicted demand
     1    2    3    4    5
1    0    0    0    0    0
2    106  88   118  131  95
3    23   21   22   18   19
4    51   73   56   51   62
5    64   70   67   53   66
6    63   55   66   59   60
7    66   48   67   58   54
8    41   42   49   45   31
9    45   38   31   37   36
10   66   73   57   59   47
11   97   79   80   87   79
12   77   86   86   83   81
13   50   67   79   76   54  ;
variable Z objective function ;
 binary variable  y(t,k),yt(t)  ;
  positive variable
                  Dt(i,t)  Order
                   s(i,t)   Shortage
                    In(i,t)  Inventory
transportation   transportation cost
Inventory        inventory cost
Shortage         shortage cost ;

Equation
eq_obj,eq24,eq25,eq26,eq27,eq28,eq291,eq292,eq301,eq302,eq311,eq312
e1
e2
e3

;

e1.. transportation=e=  sum(t,C*yt(t))     ;
e2.. Inventory =e=  sum((i,t)$(ord(i) ne 1),Ci*In(i,t)) ;
e3..  Shortage =e=  sum((i,t)$(ord(i) ne 1),Cs*d1(i)*S(i,t))  ;


eq_obj..  Z =e= transportation + Inventory + Shortage    ;


eq24(t,k,l)$(ord(t)>1 and ord(k)>ord(t) and ord(l)<ord(t))..  y(t,t)+2*y(k,l)=l=2 ;
eq25(t)..                                                     y(t,t)=e=yt(t) ;
eq26(k,t)$(ord(k)>ord(t))..                                   y(k,t)=l=yt(t) ;
eq27(t)..                                                     sum(k$(ord(t) >= ord(k)),y(t,k)) =e= 1 ;
eq28(i,t)$ (ord(i) > 1 )..                                    sum(k$(ord(t) <= ord(k)),dprime(i,k)*y(k,t)) =e= Dt(i,t) ;
eq291(i,t)$ ( ord(t)>1 and ord(i)>1)..                        In(i,t)=g= In(i,t-1)+Dt(i,t)-d(i,t)    ;
eq292(i,t)$ (ord(t)=1 and ord(i)>1)..                         In(i,t)=g= Dt(i,t)-d(i,t) ;
eq301(i,t)$(ord(t)>1and ord(i)>1)..                           S(i,t)=g= d(i,t)-Dt(i,t)-In(i,t-1) ;
eq302(i,t)$(ord(t)=1and ord(i)>1)..                           S(i,t)=g= d(i,t)-Dt(i,t) ;
eq311(i,t)$(ord(i)>1 and ord(t)>1)..                          Dt(i,t)+In(i,t-1)-d(i,t)-In(i,t)+S(i,t) =e= 0 ;
eq312(i,t)$(ord(i)>1 and ord(t)=1)..                          Dt(i,t)-d(i,t)-In(i,t)+S(i,t) =e= 0 ;

Model mode2 /all/ ;
option optcr=0
solve mode2 using MIP minimizing Z ;


