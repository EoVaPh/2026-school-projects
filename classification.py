import networkx as nx

data = []
with open("records/data.txt", "r", encoding = "utf8") as f:
    lines = f.readlines()
    for i in lines:
        pair = i.split()
        pair[2] = float(pair[2].replace("\n", ""))
        data.append(pair)

THRESHOLD = 0.5

# Создаем граф
G = nx.Graph()

# Добавляем только те ребра, которые проходят порог
valid_edges = [(u, v) for u, v, score in data if score >= THRESHOLD]
G.add_edges_from(valid_edges)
G.add_edge("8rvt", "8sbd")
# Добавляем оставшиеся одиночные узлы, которые вообще не прошли порог
all_nodes = set(u for triple in data for u in triple[:2])
G.add_nodes_from(all_nodes)

# Получаем списки семейств
families = list(nx.connected_components(G))
family = ""
with open ("records/families.txt", "w", encoding="utf8") as f:
    for i in families:
        for j in i:
            family = family + j + " "
        f.write(family)
        f.write("\n")
        family = ""