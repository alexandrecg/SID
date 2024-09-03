# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 16:09:38 2020

@author: A. Goulart
"""
value="teste123"

number = 100

class read:
    def __init__(self, val):
        print("object created with name ",val)
    
    def number(self,val):
        self.num = val
        return self.num
        
    def disp_num(self):
        print(self.num)
    
objeto = read(value)

objeto.number(number)

objeto.disp_num()

class method:
    def __init__(self,nome,valor,numero):
        self.name = nome
        self.val = valor
        
        
        B = read(nome)
        
        self.num = B.number(numero)
        
        B.disp_num
        
        
        
        print("object created with name: ", self.name, "(num: ", self.num,";val: ", self.val ,")")
        
    def __str__(self):
        return self.name
        

objeto2 = method("opa", value, number)

print(objeto2.__class__.__name__)

print(objeto2)