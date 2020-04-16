print('Oleracea')
core_oler = set([x.rstrip() for x in open('Oleracea_Pangenome_Core_GeneNames.txt')])

print('Core count:%s'%(len(core_oler)))
var_oler = set([x.rstrip() for x in open('Oleracea_Pangenome_Variable_GeneNames.txt')])
print('Variable count: %s'%(len(var_oler)))

core_in_links = set()
var_in_links = set()
for line in open('Oleracea_In_Links.txt'):
    ll = line.split()
    if ll[0] in core_oler:
        core_in_links.add(ll[0])
    if ll[0] in var_oler:
        var_in_links.add(ll[0])
print('Core in links: %s'%(len(core_in_links)))
print('Variable in links: %s'%(len(var_in_links)))

print('---------------------')
print('Rapa')
core_rapa = set([x.rstrip() for x in open('Rapa_Pangenome_WithFPSc_Core_GeneNames.txt')])

print('Core count:%s'%(len(core_rapa)))
var_rapa = set([x.rstrip() for x in open('Rapa_Pangenome_WithFPSc_Variable_GeneNames.txt')])
print('Variable count: %s'%(len(var_rapa)))

core_in_links = set()
var_in_links = set()
for line in open('Rapa_In_Links.txt'):
    ll = line.split()
    if ll[0] in core_rapa:
        core_in_links.add(ll[0])
    if ll[0] in var_rapa:
        var_in_links.add(ll[0])
print('Core in links: %s'%(len(core_in_links)))
print('Variable in links: %s'%(len(var_in_links)))

print('---------------------')

core_without_synths = set([x.rstrip() for x in open('Napus_Pangenome_Without_Synths_Core_GeneNames.txt')])
var_without_synths = set([x.rstrip() for x in open('Napus_Pangenome_Without_Synths_Variable_GeneNames.txt')])

print('Core without synths count: %s'%len(core_without_synths))
core_a = 0
core_c = 0
for a in core_without_synths:
    if 'chrA' in a:
        core_a += 1
    if 'chrC' in a:
        core_c += 1
print('Core A without synths count: %s'%core_a)
print('Core C without synths count: %s'%core_c)

core_in_links = set()
var_in_links = set()
for line in open('Napus_In_Links.txt'):
    ll = line.split()
    if ll[0] in core_without_synths:
        core_in_links.add(ll[0])
    if ll[0] in var_without_synths:
        var_in_links.add(ll[0])

print('Core in links: %s'%(len(core_in_links)))

print('Var without synths count: %s'%(len(var_without_synths)))
var_a = 0
var_c = 0
for a in var_without_synths:
    if 'chrA' in a:
        var_a += 1
    if 'chrC' in a:
        var_c += 1
print('Variable A without synths count: %s'%var_a)
print('Variable C without synths count: %s'%var_c)

print('Variable in links: %s'%len(var_in_links))

a_core_in_links = 0
c_core_in_links = 0
total_a_c = 0
for a in core_in_links:
    if 'chrA' in a or 'chrC' in a:
        total_a_c += 1
    if 'chrA' in a:
        a_core_in_links += 1
    if 'chrC' in a:
        c_core_in_links += 1 

print('Core AC in links: %s'%(total_a_c))
print('Core only A in links: %s'%a_core_in_links)
print('Core only C in links: %s'%c_core_in_links)

a_var_in_links = 0
c_var_in_links = 0
total_a_c = 0
for a in var_in_links:
    if 'chrA' in a or 'chrC' in a:
        total_a_c += 1
    if 'chrA' in a:
        a_var_in_links += 1
    if 'chrC' in a:
        c_var_in_links += 1 

print('var AC in links: %s'%(total_a_c))
print('var only A in links: %s'%a_var_in_links)
print('var only C in links: %s'%c_var_in_links)

print('--------------')
print('Now the same WITH synths')

core_with_synths = set([x.rstrip() for x in open('Napus_Pangenome_with_Synths_Core_GeneNames.txt')])
var_with_synths = set([x.rstrip() for x in open('Napus_Pangenome_with_Synths_Variable_GeneNames.txt')])

print('Core with synths count: %s'%len(core_with_synths))
core_a = 0
core_c = 0
for a in core_with_synths:
    if 'chrA' in a:
        core_a += 1
    if 'chrC' in a:
        core_c += 1
print('Core A with synths count: %s'%core_a)
print('Core C with synths count: %s'%core_c)



core_in_links = set()
var_in_links = set()
for line in open('Napus_In_Links.txt'):
    ll = line.split()
    if ll[0] in core_with_synths:
        core_in_links.add(ll[0])
    if ll[0] in var_with_synths:
        var_in_links.add(ll[0])

print('Core in links: %s'%(len(core_in_links)))

print('Var with synths count: %s'%(len(var_with_synths)))
var_a = 0
var_c = 0
for a in var_with_synths:
    if 'chrA' in a:
        var_a += 1
    if 'chrC' in a:
        var_c += 1
print('Variable A with synths count: %s'%var_a)
print('Variable C with synths count: %s'%var_c)


print('Variable in links: %s'%len(var_in_links))

a_core_in_links = 0
c_core_in_links = 0
total_a_c = 0
for a in core_in_links:
    if 'chrA' in a or 'chrC' in a:
        total_a_c += 1
    if 'chrA' in a:
        a_core_in_links += 1
    if 'chrC' in a:
        c_core_in_links += 1 

print('Core AC in links: %s'%(total_a_c))
print('Core only A in links: %s'%a_core_in_links)
print('Core only C in links: %s'%c_core_in_links)

a_var_in_links = 0
c_var_in_links = 0
total_a_c = 0
for a in var_in_links:
    if 'chrA' in a or 'chrC' in a:
        total_a_c += 1
    if 'chrA' in a:
        a_var_in_links += 1
    if 'chrC' in a:
        c_var_in_links += 1 

print('var AC in links: %s'%(total_a_c))
print('var only A in links: %s'%a_var_in_links)
print('var only C in links: %s'%c_var_in_links)
