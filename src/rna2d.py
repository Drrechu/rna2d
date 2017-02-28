# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.
from __future__ import print_function
from builtins import input
#import matplotlib


#import networkx as nx
import re
import sys
__author__ = "Dawid Rech"
__date__ = "$2017-02-06 20:25:52$"
#class  rna2D():
def rna_structure(dot_bracket_rna_file):
    """Function creates string wich associate nucleotydes in RNA sequence with 
    certein structure. Builds a graph of connections between nucleotydes and 
    writes resoults to a text file.
    
    Args:
    file: file containig RNA sequence, secondary structure in dot-bracket format 
    and RNA sequence name with ">" prefix
    
    """
    
    name_seq_dot = file_loader(dot_bracket_rna_file) 
                
    name = name_seq_dot[0]
    rna_seq = name_seq_dot[1]
    rna_dotbracket = name_seq_dot[2]
    rna_dotbracket_pseudo = name_seq_dot[3]
    connections = connection_dict_creator(rna_dotbracket)        
     
    final_seq = [] 
    
    pseudo = pseudo_handler(rna_dotbracket_pseudo)    
    duplex = brackets_handler(rna_dotbracket)
    structure_types = dots_handler(rna_dotbracket,connections)

    singleChain = structure_types[0]
    hairpins = structure_types[1]
    loop = structure_types[2]
    bulge = structure_types[3]
    junction = structure_types[4]
    
                   
    for s in singleChain: final_seq.append(s)               
    for h in hairpins: final_seq.append(h)             
    for l in loop: final_seq.append(l)             
    for b in bulge: final_seq.append(b)               
    for j in junction: final_seq.append(j)       
    for d in duplex: final_seq.append(d)
    
    final_seq.sort()
    stri = ""
    for x in final_seq:
        if x in duplex: stri += "d"*(x[1]-x[0]+1)                              
        if x in singleChain: stri += "s"*(x[1]-x[0]+1)                       
        if x in hairpins: stri += "h"*(x[1]-x[0]+1)                     
        if x in loop: stri += "l"*(x[1]-x[0]+1)                     
        if x in bulge: stri += "b"*(x[1]-x[0]+1)                      
        if x in junction: stri += "j"*(x[1]-x[0]+1) 
    
    for c in rna_dotbracket_pseudo:
        if c == "-":
            i = rna_dotbracket_pseudo.index(c)
            stri = stri[:i]+"u"+stri[i:]
    stri_pseudo=""    
    for c in pseudo:
        if pseudo[c]:
            stri_pseudo += stri[c].upper()            
        else:
            stri_pseudo += stri[c]    
    
    print (
    "====================================================================\n\
Letters in structure stands for: s-single chain, h-hairpins, l-loop, \
b-bulge, j-junction, d-duplex.\nCapital letters indicates pseudo-knots. \
Letter \"u\" represents missing structure information\n\
=====================================================================")   
    name = name.strip(">")
    print ("Name:        ",name)
    print ("RNA:         ", rna_seq)
    print ("Dot-bracket: ", rna_dotbracket_pseudo)    
    print ("Structure:   ", stri_pseudo)
    graph = graph_maker(connections,pseudo)    
    graph = str(graph)
    newLine = name+'\t'+rna_seq+'\t'+rna_dotbracket_pseudo+'\t'+stri_pseudo+'\t'+graph+'\n'
    db_append(newLine)
    
    
def db_append(new_line):
    """
    Function adds new line to a db.txt file
    """
    db=open("db.txt",'r')
    for line in db:
        if line == new_line:
            return
    db=open("db.txt",'a')    
    db.write(new_line)    
    db.close
    return
    
def file_loader(file):   
    """Function searches the file for name of the RNA sequance, sequence intself 
    and it's secondary structure in dot-bracket format. Founded strings are 
    returned along with RNA secondary structure in dot-bracket format without 
    pseudo-knots.
    
    Args:
    file: file containig RNA sequence, secondary structure in dot-bracket format 
    and RNA sequence name with ">" prefix
          
    Returns:
    list: RNA name, RNA seqence, RNA structure in dot-bracket format with pseudo-knots, 
    RNA structure in dot bracket format without pseudo-knots
    """
    f = open(file)
    name=""
    rna_dotbracket_pseudo=""
    rna_seq=""
    for line in f.readlines():
        if re.match("[)(.{}<>-]",line): rna_dotbracket_pseudo = line.strip()            
        if line[0]==">": name = line.strip()
        if re.match("[AGCUuacg]",line): rna_seq = line.strip()
    if name==""  or rna_seq=="" or rna_dotbracket_pseudo=="" : print ("Wrong file format")
    rna_dotbracket = rna_dotbracket_pseudo.replace("[",".").replace("]",".").strip()
    rna_dotbracket = rna_dotbracket.replace("{",".").replace("}",".")
    rna_dotbracket = rna_dotbracket.replace("<",".").replace(">",".").replace("-",".")
    f.flush()
    f.close()  
    return [name,rna_seq,rna_dotbracket,rna_dotbracket_pseudo]

def connection_dict_creator(RNA_in_dot_bracket):
    """
    Function creates a map of connections between non neibour nucleotydes.
    
    Args:
    RNA_in_dot_bracket (str): RNA secondary structure in dot-bracket notation, 
    without pseudo knots
    
    Returns:
    connectios (dict): Dictionary with indexes of all nucleotydes as keys and 
    indexes of bonded nucleotydes as values.
    """
    connections = {}
    stack = []
    for i, j in enumerate(RNA_in_dot_bracket):
            connections[i] = None
    for i, j in enumerate(RNA_in_dot_bracket):
        if j == '(':
            stack.append(i)
        if j == ')':
            if len(stack):
                a = stack.pop()
                b = i
                connections[a] = b
                connections[b] = a 
            else:
                print("Too many closing brackets").strip("()")
                return
        
    if len(stack):
        print ("Too many opening brackets")
        return
    return connections          
def graph_maker(Connections_Dict,Pseudo_Connections):
    """
    Function creates a graph with nucleotydes as nodes and connections between 
    them as edges
    
    Args:
    Connections_Dict (dict): Dictionary with indexes of all nucleotydes as keys and 
    indexes of bonded nucleotydes as values.
    Pseudo_Connections (dict): Dictionary with indexes of all nucleotydes as keys and 
    indexes of nucleotydes that creates pseudoknots with them as values.
    
    Returns:
    graph (dict): Dictionary with indexes of all nucleotydes as keys and lists 
    of indexes of all nucleotydes connected to them.
    """
    graph={}
    for i in range(len(Connections_Dict)):
        graph[i]=[]
        if Connections_Dict[i]:graph[i]=[Connections_Dict[i]]
    for key in graph:
        if key-1>=0:graph[key].append(key-1)
        if key<len(graph)-1:graph[key].append(key+1)
        if Pseudo_Connections[key]:graph[key].append(Pseudo_Connections[key])

    return graph

def brackets_handler(RNA_in_dot_bracket):
    """
    Function recognizes duplex structure in RNA secondary structure
    
    Args: 
    
    RNA_in_dot_bracket (str): RNA secondary structure in dot-bracket notation, 
    without pseudo-knots
    
    Returns: duplex (list): List of tuples, every tuple has twe elements, 
    index of begining end end of duplex.    
    """
    stack = re.compile('[()]{1,}')
    duplex = []
    for s in stack.finditer(RNA_in_dot_bracket):
        duplex.append(s.span())
    for i ,d in enumerate(duplex):
        duplex[i]=(d[0],d[1]-1)        
    return duplex   
    

def dots_handler(RNA_in_dot_bracket,Connestions_Dict):
    """
    Function recognizes structures in RNA seconondary stucture. Recognized structures are:
    single chain, hairpin, loop, bulge and junction.
    
    Args:
    RNA_in_dot_bracket (str): RNA secondary structure in dot-bracket notation, 
    without pseudo knots
    
    Connections_Dict (dict): Dictionary with indexes of all nucleotydes as keys and 
    indexes of bonded nucleotydes as values.
    
    Returns:
    List of lists of tuples, every list is made of tuples, every tuple hes two 
    elements, the begining and the end of certein strucure.    
    """
    singleChain = []
    hairpins = []
    loop = []
    bulge = []
    junction = []    
    dots = []

    dot = re.compile('[.]{1,}')
    for dt in dot.finditer(RNA_in_dot_bracket):
        dots.append(dt.span()) 
   
    for dt in dots:
        start = dt[0] 
        stop = dt[1]        
        if start == 0 or stop>=len(RNA_in_dot_bracket)-1:
            singleChain.append((start,stop-1))
            continue
        if (RNA_in_dot_bracket[start-1]==')' and RNA_in_dot_bracket[stop]=='('):
            junction.append((start, stop-1))
            continue
       
        if (RNA_in_dot_bracket[start-1]=='(' and RNA_in_dot_bracket[stop]==')') :
           hairpins.append((start, stop-1))
           continue
           
        if (RNA_in_dot_bracket[start-1]=='(' and RNA_in_dot_bracket[stop]=='(') :
           
            if Connestions_Dict.get(start-1,  'None') -  Connestions_Dict.get(stop,  'None') == 1:
                bulge.append((start, stop-1))
            else:
                loop.append((start, stop-1))
          
            continue
        if (RNA_in_dot_bracket[start-1]==')' and RNA_in_dot_bracket[stop]==')'):
            if Connestions_Dict.get(start-1,  'None') -  Connestions_Dict.get(stop,  'None') == 1:
                bulge.append((start, stop-1))
            else:
                loop.append((start, stop-1))
            continue
            
    for cr in junction:
        start = cr[0]        
        for pet in loop:
            if pet[1] == Connestions_Dict.get(start-1,  'None') - 1: 
                junction.append((pet[0],pet[1]))
                loop.remove((pet[0],pet[1]))
    junction.sort  
    return [singleChain,hairpins,loop,bulge,junction]

def pseudo_handler(RNA_in_dot_bracket):
    """Function creates dictionary with information about pseudo-knots.
    Args:
    RNA_in_dot_bracket (str): RNA secondary structure in dot-bracket notation, 
    with pseudo knots.
    
    Returns:
    pseudo (dict): Dictionary with indexes of all nucleotydes as keys and 
    indexes of bonded nucleotydes creating pseudo-knots as values.
    """
    
    pseudo = {}
    stack = []
    stack2 = []
    stack3 = []
    for i, j in enumerate(RNA_in_dot_bracket):
            pseudo[i] = None
    for i, j in enumerate(RNA_in_dot_bracket):
        if j == '[':stack.append(i)
        if j == '{':stack2.append(i)
        if j == '<':stack3.append(i)
            
        if j == ']':
            if len(stack):
                a = stack.pop()
                b = i
                pseudo[a] = b
                pseudo[b] = a
            else:
                print("Too many closing brackets")
                return
        if j == '}':
            if len(stack2):
                a = stack2.pop()
                b = i
                pseudo[a] = b
                pseudo[b] = a
            else:
                print("Too many closing brackets")
                return
        if j == '>':
            if len(stack3):
                a = stack3.pop()
                b = i
                pseudo[a] = b
                pseudo[b] = a
            else:
                print("Too many closing brackets")
                return        
    if len(stack) or len(stack2) or len(stack3):
        print("Too many opening brackets")
        return
    return pseudo          
if __name__ == "__main__":
    try:
        rna_structure(sys.argv[1])
    except  IndexError:      
        file = input("Please type path to a file you want to processed. \t")       
        rna_structure(file)    
    input("Enter to terminate. Output has been added to the db.txt file")
    sys.exit()
   
     
