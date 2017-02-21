# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

__author__ = "User"
__date__ = "$2017-02-06 20:25:52$"
import re

def rna_structure(dot_bracket_rna_file):
    """Main Function"""
    
    name_seq_dot = file_loader(dot_bracket_rna_file) 
                
    name = name_seq_dot[0]
    rna_seq = name_seq_dot[1]
    rna_dotbracket = name_seq_dot[2]
    print rna_dotbracket
    connections = connection_dict_creator(rna_dotbracket)        
     
    final_seq = [] 
        
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
                                         
    pseudo = pseudo_handler(rna_dotbracket)
    print "speudo", pseudo
    print "Name: ",name.strip(">")
    print "RNA:  ", rna_seq
    print "db:   ", rna_dotbracket    
    print "stri: ", stri
    graph_maker(connections,rna_dotbracket)
    
    f=open("db.txt",'w') 
    f.write(name + rna_seq +rna_dotbracket + stri)
    f.closed
    
def file_loader(file):   
    """File parser"""
    f = open(file)
    rna_dotbracket_pseudo=""
    for line in f.readlines():
        if re.match("[)(.]",line): rna_dotbracket_pseudo += line            
        if line[0]==">": name = line.strip()
        if re.match("[AGCU]",line): rna_seq = line.strip()
    if name=="" or rna_seq=="" or rna_dotbracket_pseudo=="" : print "Zly format pliku"
    rna_dotbracket = rna_dotbracket_pseudo.replace("[",".").replace("]",".")
    f.closed 
    
    return [name,rna_seq,rna_dotbracket,rna_dotbracket_pseudo]

def connection_dict_creator(RNA_in_dot_bracket):
    """Function creates a map of connections between non neibour nucleotydes"""
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
                print("Too many closing brackets")
                return
        
    if len(stack):
        print("Too many opening brackets")
        return
    
    return connections          
def graph_maker(Connections_Dict,RNA_in_dot_bracket):
    for i in range(len(Connections_Dict)):
        Connections_Dict[i]=[(Connections_Dict[i])]
    graph=Connections_Dict
    for key in graph:
        if key-1>=0:graph[key].append(key-1)
        if key<len(graph)-1:graph[key].append(key+1)
    print graph
    return graph

def brackets_handler(RNA_in_dot_bracket):
    stack = re.compile('[()]{1,}')
    duplex = []
    #new_duplex= []
    for s in stack.finditer(RNA_in_dot_bracket):
        duplex.append(s.span())
    for i ,d in enumerate(duplex):
        duplex[i]=(d[0],d[1]-1)        
    return duplex   
    

def dots_handler(RNA_in_dot_bracket,Connestions_Dict):
    
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
        start = dt[0] #pierwszy z pary kropek
        stop = dt[1] #drugi z pary kropek           
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
            if pet[1] == Connestions_Dict.get(start-1,  'None') - 1: #koniec petli to skrzyzowanie
                junction.append((pet[0],pet[1]))
                loop.remove((pet[0],pet[1]))
    junction.sort  
    
    return [singleChain,hairpins,loop,bulge,junction,loop]

def pseudo_handler(RNA_in_dot_bracket):
    """Function returns pseudo"""
    
    pseudo = {}
    stack = []
    for i, j in enumerate(RNA_in_dot_bracket):
            pseudo[i] = None
    for i, j in enumerate(RNA_in_dot_bracket):
        if j == '[' or j == '{' :
            stack.append(i)
        if j == ']' or j == '}':
            if len(stack):
                a = stack.pop()
                b = i
                pseudo[a] = b
                pseudo[b] = a
            else:
                print("Too many closing brackets")
                return
        
    if len(stack):
        print("Too many opening brackets")
        return
    
    return pseudo          
if __name__ == "__main__":
         
     #elementyStruktury('..((((..((..))..))((..((..))..))))...')
     #elementyStruktury('..(((.((((((...)))....)))....)))')
     #elementyStruktury('..((((((...)))(((...)))(((...)))(((...)))(((...)))(((...)))(((...)))))).')
     #elementyStruktury('..((((((...)))(((...)))..(((...)))(((...)))(((...)))(((...)))(((...)))))).')
     #elementyStruktury('..((((((...)))(((..((..(())..))..((..(())..))..)))..(((...)))(((...)))(((...)))(((...)))(((...)))))).')
     #elementyStruktury('..(((((()))(((...)))(((...)))))).')
     #elementyStruktury('..((..(((..((..))..)))..))..')
     #elementyStruktury('..((..{(((((..))..))))..))..')
     #elementyStruktury('..((..(((..((..)))))..))..')
     #elementyStruktury('(((..(((...)))...(((...)))..)))')
     
     #elementyStruktury('((..((..((..((..((..((..((..))..((..((..((..))..))..))..))..))..))..((..))..((..((..((..((..((..))..))..((..((..))..))..))..))..))..((..((..((..((..))..((..))..))..))..))..))..))..))')
     
     #elementyStruktury('((((((((..((((....((..((((..((..(((((.....)))))....(((((..((..((((....))))..))....)))))..)).....))..))..))....(((((....)))..))..((...((..((..((((((....((((....))))..))))))..((..((....))..))....))..))..))..((..((..((((..((((...))))..((((...))))..))))..))..))..))..))..))))))))......')  
     #elementyStruktury('..(((...(((.....))))))...')
     #elementyStruktury('((..((..((..))..((..))..))..((..))..))')
     #elementyStruktury('((..((..((..))((..))..))..((..))..))')
     #elementyStruktury('..((((..(((((.(((((((((....)))))))))..)))))....((((((((....((((.(((((....))))).)))).)))))))).))))')
   
     rna_structure('file.txt')

     #PDB_00176
     #rna_structure('.(((((...((((.)))).))))).(((((...((((.)))).))))).(((((...((((.)))).))))).(((((...((((.)))).)))))')
     #
     
   
     