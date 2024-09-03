# SID
**SID** stands for **Swirl Injector Design**. This project originated from some of the codes I've developed for my undergraduate thesis.
 
![alt text](https://github.com/alexandrecg/SID/blob/main/RD0110_homepage.png "RD0110 Main Injector Initial Modeling")

![alt text](https://github.com/alexandrecg/SID/blob/main/RD0110_F_mass_flow%20-%20homepage.png "RD0110 Fuel Injector Mass Flow Modeling")

The original concept of the software developed in this project contemplates three tools.

**1)** Design from the input of basic injector properties:
  You can use this information to get a corresponding injector geometry using basic engine parameters and a proposed injection plate design.

**2)** Design from the input of injector geometry:
  Inputing injector geometry, flow and spray characteristics are predicted.
  
**3)** Manufacturing tolerances evaluation:
  When injector geometry and manufacturing tolerances are input, an injector operation envelope is predicted.
  
  
Tool 1) is capable of presenting monopropellant and bipropellant designs, following the method presented by Bazarov in reference [1]. However, in its development and study, other methods, like the one presented by Bayvel [2], are also compared, especially concerning a difference in friction factors. Both of them have at their core the works of Abramovich and Kliachko, and the applicability of these theories is well discussed in the work of Khavkin [3].

To guarantee an adequate degree of correspondence with the data calculated by the proposed methods, Bipropellant injectors are always checked for hydraulic independence between their stages; furthermore, internal impingement is not contemplated for the same reason.



**References:**

[1]: Vladimir Bazarov, Vigor Yang, Puneesh Puri (2004) 'Design and Dynamics of Jet and Swirl Injectors', in Vigor Yang, Mohammed Habiballah, James Hulka, Michael Popp, Paul Zarchan (ed.) Liquid Rocket Thrust Chambers: Aspects of Modeling, Analysis, and Design. : American Institute of Aeronautics and Astronautics, pp. 19-103.

[2]: L. Bayvel, Z. Orzechowski (1993) Liquid Atomization, Taylor & Francis Books.

[3]: Yuriy I. Khavkin (2004) Theory and Practice of Swirl Atomizers, Taylor & Francis Books.

[--]: A.  Alves,  “Estudo  e  Desenvolvimento  de  um  Sistema  de  Injeção  Centrífugo Bipropelente  Utilizado  e  Motor  Foguete  a  Propelente  Líquido.,”  Instituto Tecnológico de Aeronáutica, 2008.

[--]: A. H. Lefebvre e V. G. Mcdonell., Atomization and Sprays. 2th ed., CRC Press, Taylor & Francis Group, 2017. 
