import networkx as nx
import sys
import heapq
from preprocessing import get_reads


class DeBruijn:
    def __init__(self, reads):
        self.reads = reads
        self.k = len(reads[0])
        self.graph = None
        self.dot_file = None

    def construct_graph(self):
        self.graph = nx.DiGraph()

        for ele in self.reads:
            if not self.graph.has_edge(ele[:self.k-1], ele[1:]):
                self.graph.add_edge(ele[:self.k-1], ele[1:], content=ele[-1] , cnt=1)
            else:
                self.graph.edges[ele[:self.k-1], ele[1:]]["cnt"] += 1
        labels = set()
        nx.set_node_attributes(self.graph, labels, "labels")
        # print(len(self.graph.nodes))
        # print(len(self.graph.edges))

    def merge_singleton(self):
        '''
        Singleton:
            A --> B
            B only one in degree
            A only one out degree
        '''
        # update_content = self.graph["TATA"]["ATAC"]["content"]

        # self.graph = nx.contracted_nodes(self.graph, "TATA", "ATAC", self_loops = False)
        # self.graph = nx.relabel_nodes(self.graph, {"TATA": "TATA"+update_content})
        # print(list(self.graph.predecessors))
        # self.graph["CTAT"]["TATAC"]["content"] = "test"

        # print(nx.get_edge_attributes(self.graph, "content")[("CTAT","TATAC")])
        cnt = 0
        # print(list(self.graph.predecessors("CGCT")))
        # while cnt <= len(self.graph.nodes):
        init = list(nx.dfs_preorder_nodes(self.graph))
        # print(init)
        for i, node in enumerate(init):
            # print(f'--{i}--')
            # print("node:", node, self.graph.nodes())
            if node not in self.graph:
                continue
            # print('Get predecessor')
            prev_nodes = list(self.graph.predecessors(node))
            # print("f", node, prev_nodes)
            if len(prev_nodes) == 1 and self.graph.out_degree(prev_nodes[0]) == 1:
                # print('Find Singleton ...')
                prev_node = prev_nodes[0]
                update_content = self.graph[prev_node][node]["content"]
                # print(f'Graph Before Merge: Node: {len(self.graph.nodes())}, {len(self.graph.edges())}')
                self.graph = nx.contracted_nodes(self.graph, prev_node, node, self_loops = False)
                # print(f'Graph After Merge: Node: {len(self.graph.nodes())}, {len(self.graph.edges())}')
                # print('Finish Merge ...')

                new_prev_node = prev_node+update_content
                # print(new_prev_node)
                # print(prev_node)
                self.graph = nx.relabel_nodes(self.graph, {prev_node: new_prev_node})
                # print('Finish Relable ...')

                # update prev node in degree content
                for pre_prev_node in list(self.graph.predecessors(new_prev_node)):
                    self.graph[pre_prev_node][new_prev_node]["content"] += update_content
                # print('Finish Update ...')
                # print("check", prev_node, node)
                # cnt = 0
                # break
            # else:
            #     cnt += 1
        # print(self.graph.edges)
        # nx.draw(self.graph, with_labels=True)
        # plt.show()
    def detect_cycle(self):
        print("detect cycle", nx.find_cycle(self.graph))

    def get_start_node(self):
        nodes = self.graph.nodes
        start_nodes = []
        for node in nodes:
            if self.graph.in_degree(node) == 0:
                start_nodes.append(node)
        return start_nodes


    def find_path(self, start_nodes):
        res = []
        self.imax = 0

        def cycle_dfs(root, path, seen):
            if self.graph.out_degree(root) == 0 or root in seen:
                self.imax = max(self.imax, len(path))
                heapq.heappush(res, (len(path), path))
                # res.append(path)
                return res

            next_nodes = list(self.graph.successors(root))
            seen.add(root)
            for next in next_nodes:
                content = self.graph[root][next]["content"]
                # print(root, next, content)
                cycle_dfs(next, path+content, seen)

        def dfs(root, path):
            if self.graph.out_degree(root) == 0:
                self.imax = max(self.imax, len(path))
                res.append(path)
                return res
            next_nodes = list(self.graph.successors(root))

            for next in next_nodes:
                content = self.graph[root][next]["content"]
                # print(root, next, content)
                dfs(next, path+content)

        def stack_dfs(root, path):
            WHITE, GRAY = 0, 1
            stack = [(WHITE, root, path)]

            while stack:
                # print("stack", len(stack))
                color, node, cur_path = stack.pop()
                # print("fuck", node, cur_path)
                if self.graph.out_degree(node) == 0:
                    self.imax = max(self.imax, len(cur_path))
                    res.append(cur_path)
                    # print(res)
                    continue
                if color == WHITE:
                    next_nodes = list(self.graph.successors(node))
                    # print(next_nodes, self.graph.nodes)

                    for next in next_nodes:
                        content = self.graph[node][next]["content"]
                        stack.append((WHITE, next, cur_path+content))
                    stack.append((GRAY, node, cur_path))
            # print(stack)
            return res




        for idx, root in enumerate(start_nodes):
            # dfs(root, root)
            cycle_dfs(root, root, set())

            # res.append(path)
        update_res = []
        for length, contig in res:
            if length == self.imax:
                update_res.append(contig)
        return sorted(update_res), res

    def print_out(self, res):
        for ele in res:
            print(ele)

        
    def formulate_dot(self, read_file):
        edges = self.graph.edges

        with open("debruijn.dot", "w") as f:
            f.write("digraph {\n")
            for iso in list(nx.isolates(self.graph)):
                f.write("{0};\n".format(iso))
            for node1, node2 in edges:
                edge = node1 + self.graph[node1][node2]["content"]
                f.write('{0} -> {1} [label="{2}"];\n'.format(node1, node2, edge))
            f.write("}")



def main(argv):
    read_file = argv[1]
    reads = get_reads(read_file)
    graph = DeBruijn(reads)
    graph.construct_graph()
    graph.merge_singleton()
    start_nodes = graph.get_start_node()
    res, _ = graph.find_path(start_nodes)
    graph.print_out(res)
    graph.formulate_dot(argv[1])






if __name__ == "__main__":
    main(sys.argv)

