# SID
SID stands for Swirl Injector Design. This project originated from some of the codes I've developed for my undergraduate thesis.


The original concept of the software developed in this project contemplates 3 tools.

1) Design from input of basic injector properties:
  In general, if you have the basic engine parameters and a proposed injection plate design, you can use this information to get a corresponding injector geometry.

2) Design from input of injector geometry:
  Inputing injector geometry, flow and spray characteristics are predicted.
  
3) Manufacturing tolerances evaluation:
  Inputing injector geometry and manufacturing tolerances, a injector operation envelope is predicted.
  
  
Tool 1) is capable of presenting monopropellant and bipropellant designs, following the method presented by Bazarov in reference [1], but in it's development and study, other methods like the one presented by Bayvel [2] is also compared, specially concerning a difference in friction factors. Both of them have in it's core the works of Abramovich and Kliachko and applicability of these theories are well discussed in the work of Khavkin [3].

To garantee an adequate degree of correspondece on the data calculated by the proposed methods, Bipropellant injectors are always checked for hydraulic independence between its stages, furthermore, internal impingement is not contemplated for the same reason.



References:

[1]: Vladimir Bazarov, Vigor Yang, Puneesh Puri (2004) 'Design and Dynamics of Jet and Swirl Injectors', in Vigor Yang, Mohammed Habiballah, James Hulka, Michael Popp, Paul Zarchan (ed.) Liquid Rocket Thrust Chambers:Aspects of Modeling, Analysis,and Design. : American Institute of Aeronautics and Astronautics, pp. 19-103.

[2]: L. Bayvel, Z. Orzechowski (1993) Liquid Atomization, : Taylor & Francis Books.

[3]: Yuriy I. Khavkin (2004) Theory and Practice of Swirl Atomizers, : Taylor & Francis Books.

[-]: A.  Alves,  “Estudo  e  Desenvolvimento  de  um  Sistema  de  Injeção  Centrífugo Bipropelente  Utilizado  e  Motor  Foguete  a  Propelente  Líquido.,”  Instituto Tecnológico de Aeronáutica, 2008.

[-]: A. H. Lefebvre e V. G. Mcdonell., Atomization and Sprays. 2th ed., CRC Press, Taylor & Francis Group, 2017. 
