import networkx as nx


def fonction(a, b, **option):
    coucou = option.get('coucou', "default coucou")
    coucou2 = option.get('coucou2', "default coucou2")
    print(coucou)
    print(coucou2)
    return a+b


def print_my_list(title_list, title_dict, *my_list, **my_dict):
    # Dans ce namespace, my_dict est un dictionnaire et my_list une liste

    print('ici les elements de la liste %s' % (title_list))
    print(my_list)
    for idx, value in enumerate(my_list):
        print('Index %02d: %s' % (idx, value))

    print('voici les elements du dictionnaire %s' % (title_dict))
    for key, value in my_dict.items():
        print('Element %s: %s' % (key, value))


print_my_list('liste 0', 'dict 0', 0,
              zero=0, un=1, deux=2, trois=3, quatre=4, cinq=5)


print(fonction(5, 3, coucou="coucou"))
args = {"coucou": "yo", "coucou2": "yo2"}
print(fonction(5, 3, **args))

G = nx.Graph()
G.add_edge(0, 1)
print(type(G))

print((5,))
