import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.*;

import static java.lang.Math.pow;

public class Calculator {
    private static final HashSet<String> buff = new HashSet<>(); // 缓冲区
    static HashMap<String, Double> netHashMap; // 原始网络数据 对应net.txt文件 存储为gene1:gene2 value
    static ArrayList<String> geneList; // 基因列表 本应来自net.txt文件，此处由于进一步筛选因此直接从geneList.txt读入
    static ArrayList<ArrayList<Double>> netMatrix; // 原始网络矩阵,数值来自netHashMap，行列坐标来自geneList
    static HashMap<String, Integer> geneIndexMap; // 基因索引映射表
    static double[][] geneSimilarityMatrix; // 随机游走的基因相关性得分矩阵
    static HashMap<String, ArrayList<String>> sonsBuff; // 直接后代子节点缓冲区
    static ArrayList<Annotation> annotationList; // gene.gaf 注释文件数据
    static HashMap<String, ArrayList<String>> termToGeneList; // 本体术语-->术语注释的基因集间的映射表
    static HashMap<String, ArrayList<String>> idParentsHashMap; // 术语网络拓扑图父节点 哈希表: key=本体id value=本体结构信息
    static HashMap<String, ArrayList<String>> idChildrenHashMap; // 术语网络拓扑图后代节点 哈希表: key=本体id value=本体结构信息
    static HashMap<String, Set<String>> childrenGenesMap; // 直接后代基因集合 哈希表: key=基因 value=直接后代基因集合
    static HashMap<String, ArrayList<String>> idSonHashMap; // 术语网络拓扑图直接儿子节点 哈希表: key=本体id value=本体结构信息
    static HashMap<String, Term> idTermHashMap; // 本体术语项 哈希表: key=基因本体id value=基因本体结构信息
    static HashMap<String, ArrayList<Annotation>> idAnnotationHashMap; // 注释信息数据 哈希表: key=本体ID value=注释信息
    static HashMap<String, ArrayList<String>> geneToTermList;// 基因-->注释基因的本体术语集合间的映射表
    static HashMap<String, HashSet<String>> ecNumGeneHashMap; // 基因证据编号 哈希表: key=证据编码 value=该ec编码的基因集
    static HashMap<String, Double> similarityResult; // 相似度计算结果
    private static boolean isInit = false; // 是否初始化标志
    // 基因索引映射表

    // Part1
    // 判断文件是否存在，并创建文件
    public static void createFileIfNotExists(Path path) throws IOException {
        if (!Files.exists(path)) {
            Files.createFile(path);
        }
    }

    // 根据参数选择文件读写路径
    public static ArrayList<String> getFilePaths(String gene_type) {
        ArrayList<String> filePaths = new ArrayList<>();
        if (Objects.equals(gene_type, "Yeast")) {
            filePaths.add("buf/Yeast/YeastNet.v3.txt");
            filePaths.add("buf/Yeast/YeastNet.v3.benchmark.txt");
            filePaths.add("buf/Yeast/sgd.gaf");
            filePaths.add("buf/Yeast/biochemical_pathways.tab");
        } else if (Objects.equals(gene_type, "Human")) {
            filePaths.add("buf/Human/HumanNet-XC.tsv");
            filePaths.add("buf/Human/HumanNet-GSP.tsv");
            filePaths.add("buf/Human/goa_human.gaf");
            filePaths.add("buf/Human/PubChem_pathway_text_human.csv");
        } else if (Objects.equals(gene_type, "Arabidopsis_thaliana")) {
            filePaths.add("buf/Arabidopsis_thaliana/AraNet.txt");
            filePaths.add("buf/Arabidopsis_thaliana/AraNet_GS.txt");
            filePaths.add("buf/Arabidopsis_thaliana/tair.gaf");
            filePaths.add("buf/Arabidopsis_thaliana/aracyc_pathways.20230103");
        }
        return filePaths;
    }

    // 计算并存储原始基因网络数据矩阵 netMatrix
    private static void getNetMatrix(String gene_type) throws IOException {
        // 根据参数修改文件路径
        final ArrayList<String> filePaths = getFilePaths(gene_type);
        // 初始化 网络数据哈希表
        netHashMap = new HashMap<>();
        for (FunctionNet edge : Reader.readNetFile(filePaths.get(0))) {
            netHashMap.put(edge.gene1 + ":" + edge.gene2, edge.value / 10.0D); //做归一化处理
        }
        // 读取 黄金基因列表
        HashSet<String> goldGene = new HashSet<>();
        for (String line : Files.readAllLines(Paths.get(filePaths.get(1)))) {
            goldGene.add(line.split("\t")[0]);
            goldGene.add(line.split("\t")[1]);
        }
        geneList = new ArrayList<>(goldGene);
        // 存储 基因列表
        StringBuilder sg = new StringBuilder();
        for (String gene : geneList) {
            sg.append(gene).append("\n");
        }
        final Path genePath = Paths.get("temp/" + gene_type + "/geneList.txt");
        createFileIfNotExists(genePath);
        Files.write(genePath, sg.toString().getBytes());

        // 计算 原始网络矩阵
        netMatrix = new ArrayList<>();
        for (int k = 0; k < geneList.size(); k++) {
            netMatrix.add(new ArrayList<>());
            netMatrix.get(k).ensureCapacity(geneList.size());
            for (String s : geneList) {
                String key = geneList.get(k) + ":" + s;
                if (netHashMap.containsKey(key)) {
                    netMatrix.get(k).add(netHashMap.get(key));
                }
            }
        }
        // 存储 netMatrix矩阵
        StringBuilder sb = new StringBuilder();
        for (ArrayList<Double> line : netMatrix) {
            for (double value : line) {
                sb.append(value).append("\t");
            }
            sb.append("\n");
        }
        final Path netMatrixPath = Paths.get("temp/" + gene_type + "/netMatrix.txt");
        createFileIfNotExists(netMatrixPath);
        Files.write(netMatrixPath, sb.toString().getBytes());
    }

    // 计算权重矩阵
    private static ArrayList<ArrayList<Double>> getWMatrix() throws Exception {
        double n = netMatrix.size();
        double sita_2 = parameterVariable(n);
        ArrayList<ArrayList<Double>> wM = new ArrayList<>();
        wM.ensureCapacity(geneList.size());
        for (int i = 0; i < n; i++) {
            wM.add(new ArrayList<>());
            wM.get(i).ensureCapacity(geneList.size());
            for (ArrayList<Double> matrix : netMatrix) {
                wM.get(i).add(Math.exp(-(pow(euclideanDistance(netMatrix.get(i), matrix), 2) / sita_2)));
            }
        }
        return wM;
    }

    // 计算参数变量的平方
    private static double parameterVariable(double n) throws Exception {
        double sum = 0.0D;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i != j)
                    sum += euclideanDistance(netMatrix.get(i), netMatrix.get(j));
            }
        }
        sum /= (n * (n - 1));
        sum *= sum;
        return sum;
    }

    // 计算两向量间的距离，欧氏距离
    private static double euclideanDistance(ArrayList<Double> vector1, ArrayList<Double> vector2) throws Exception {
        double value = 0.0D;
        if (vector1.size() != vector2.size())
            throw new Exception("向量不等长");
        for (int i = 0; i < vector1.size(); i++)
            value += Math.pow((vector1.get(i) - vector2.get(i)), 2);
        value = Math.sqrt(value);
        return value;
    }

    // 计算概率转移矩阵 pMatrix
    private static ArrayList<ArrayList<Double>> getPMatrix() throws Exception {
        ArrayList<ArrayList<Double>> wMatrix = getWMatrix();
        ArrayList<ArrayList<Double>> pMatrix = new ArrayList<>();
        pMatrix.ensureCapacity(geneList.size());
        for (int i = 0; i < netMatrix.size(); i++) {
            pMatrix.add(new ArrayList<>());
            pMatrix.get(i).ensureCapacity(geneList.size());
            for (int j = 0; j < netMatrix.size(); j++) {
                if (netMatrix.get(i).get(j) != 0.0D) {
                    double value = wMatrix.get(i).get(j);
                    double sum = 0.0D;
                    for (int k = 0; k < netMatrix.size(); k++)
                        sum += wMatrix.get(i).get(k);
                    value /= sum;
                    pMatrix.get(i).add(value);
                } else {
                    pMatrix.get(i).add(0.0D);
                }
            }
        }
        return pMatrix;
    }

    // 存储概率转移矩阵 pMatrix
    private static void storagePTMatrix(String gene_type) throws Exception {
        // 根据参数修改文件路径
        String filePath = "temp/" + gene_type + "/pMatrix.txt";
        ArrayList<ArrayList<Double>> pm = getPMatrix();
        StringBuilder sb = new StringBuilder();
        for (ArrayList<Double> line : pm) {
            for (double value : line) {
                sb.append(value).append("\t");
            }
            sb.append("\n");
        }
        final Path path = Paths.get(filePath);
        createFileIfNotExists(path);
        Files.write(path, sb.toString().getBytes());
    }

    // 计算随机游走矩阵即基因之间的相关性得分矩阵*****************************************
    public static void getRWRMatrix(String gene_type, double arg) throws Exception {
        getNetMatrix(gene_type);
        storagePTMatrix(gene_type);
        RandWalk.walk(gene_type, arg);
    }

    // Part2
    // 数据初始化
    // 获取后代节点
    private static void getChildren(String goID) {
        if (!sonsBuff.get(goID).isEmpty()) {
            for (String go : sonsBuff.get(goID)) {
                getChildren(go);
            }
        }
        buff.add(goID);
    }

    // 基因本体树构建
    private static void createGOTree(ArrayList<Term> obo_terms, String gene_type) throws Exception {
        // 生成 son.txt
        sonsBuff = new HashMap<>();
        for (Term term : obo_terms) {
            if (term.isObsolete) {
                continue;
            }
            sonsBuff.putIfAbsent(term.id, new ArrayList<>());
            // 添加子节点
            for (String parent : term.isID) {
                sonsBuff.putIfAbsent(parent, new ArrayList<>());
                sonsBuff.get(parent).add(term.id);
            }
            for (String parent : term.partID) {
                sonsBuff.putIfAbsent(parent, new ArrayList<>());
                sonsBuff.get(parent).add(term.id);
            }
        }
        // 写入 son.txt
        final Path sonPath = Paths.get("temp/" + gene_type + "/son.txt");
        createFileIfNotExists(sonPath);
        try (BufferedWriter writer = Files.newBufferedWriter(sonPath)) {
            WriteChildOrSon(sonsBuff, writer);
        }
        // 生成 child.txt
        HashMap<String, ArrayList<String>> childs = new HashMap<>();
        for (String parent : sonsBuff.keySet()) {
            childs.putIfAbsent(parent, new ArrayList<>());
            for (String child : sonsBuff.get(parent)) {
                getChildren(child);
            }
            childs.get(parent).addAll(buff);
            buff.clear();
        }
        // 写入 child.txt
        final Path childPath = Paths.get("temp/" + gene_type + "/child.txt");
        createFileIfNotExists(childPath);
        try (BufferedWriter writer = Files.newBufferedWriter(childPath)) {
            WriteChildOrSon(childs, writer);
        }
    }

    // 写入 child.txt 或 son.txt
    private static void WriteChildOrSon(HashMap<String, ArrayList<String>> childs, BufferedWriter writer) throws IOException {
        for (Map.Entry<String, ArrayList<String>> entry : childs.entrySet()) {
            String parent = entry.getKey();
            String children = String.join("\t", entry.getValue());
            writer.write(parent + ":" + children);
            writer.newLine();
        }
    }

    // 初始化子节点或后代节点
    private static void initChildOrSon(String line, HashMap<String, ArrayList<String>> idChildrenHashMap) {
        ArrayList<String> childs = new ArrayList<>();
        if (line.split(":").length > 1) {
            childs.addAll(Arrays.asList(line.split(":")[1].split("\t")));
        }
        idChildrenHashMap.put(line.split(":")[0], childs);
    }

    // 相关数据初始化
    public static void data_initializer(String gene_type, double alpha) throws Exception {
        if (isInit) {
            return;
        }
        // 根据参数修改文件路径
        ArrayList<String> filePaths = getFilePaths(gene_type);
        ArrayList<Term> obo_terms = Reader.readOboFile("buf/go-basic.obo");
        // 构建基因本体树相关数据
        createGOTree(obo_terms, gene_type);

        // 本体术语项 哈希表 初始化
        idTermHashMap = new HashMap<>();
        for (Term term : obo_terms) {
            if (!term.isObsolete) {
                idTermHashMap.put(term.id, term);
            }
        }

        // 术语网络拓扑图父节点 哈希表 初始化
        idParentsHashMap = new HashMap<>();
        for (String id : idTermHashMap.keySet()) {
            ArrayList<String> parents = new ArrayList<>();
            idParentsHashMap.put(id, parents);
            // 已经添加的术语
            HashSet<String> added = new HashSet<>();
            // 等待添加的术语
            ArrayList<String> adding = new ArrayList<>(idTermHashMap.get(id).partID);
            adding.addAll(idTermHashMap.get(id).isID);
            while (!adding.isEmpty()) {
                String current = adding.remove(0);  // 直接使用remove减少操作
                if (added.add(current)) {  // 如果未添加则返回true
                    parents.add(current);
                    Term t = idTermHashMap.get(current);
                    if (t != null) {  // 检查术语是否存在
                        if (t.partID != null) adding.addAll(t.partID);
                        if (t.isID != null) adding.addAll(t.isID);
                    }
                }
            }
        }

        // 术语网络拓扑图后代节点 哈希表 初始化
        idChildrenHashMap = new HashMap<>();
        for (String line : Files.readAllLines(Paths.get("temp/" + gene_type + "/child.txt"))) {
            initChildOrSon(line, idChildrenHashMap);
        }

        // 术语网络拓扑图直接儿子节点 哈希表 初始化
        idSonHashMap = new HashMap<>();
        for (String line : Files.readAllLines(Paths.get("temp/" + gene_type + "/son.txt"))) {
            initChildOrSon(line, idSonHashMap);
        }

        // 注释文件数据 gene.gaf
        annotationList = Reader.readGafFile(filePaths.get(2), gene_type);

        // ecNumGeneHashMap 基因证据编号 哈希表 初始化
        ecNumGeneHashMap = LfcPair.getLfc(gene_type);

        // 初始化注释信息数据哈希表、术语与基因集的映射表、基因与术语的映射表
        idAnnotationHashMap = new HashMap<>();
        termToGeneList = new HashMap<>();
        geneToTermList = new HashMap<>();

        for (Annotation annotation : annotationList) {
            // 填充 idAnnotationHashMap
            idAnnotationHashMap.computeIfAbsent(annotation.goID, k -> new ArrayList<>()).add(annotation);

            // 填充 termToGeneList
            termToGeneList.computeIfAbsent(annotation.goID, k -> new ArrayList<>()).add(annotation.geneID);
            if (annotation.synonym != null) {
                termToGeneList.get(annotation.goID).addAll(annotation.synonym);
            }

            // 填充 geneToTermList
            geneToTermList.computeIfAbsent(annotation.geneID, k -> new ArrayList<>()).add(annotation.goID);
            if (annotation.synonym != null) {
                for (String gene : annotation.synonym) {
                    geneToTermList.computeIfAbsent(gene, k -> new ArrayList<>()).add(annotation.goID);
                }
            }
        }

        // 去重处理
        termToGeneList.forEach((term, genes) -> arrayDistinct(genes));
        geneToTermList.forEach((gene, terms) -> arrayDistinct(terms));

        // 在初始化时填充所有子节点基因
        childrenGenesMap = new HashMap<>();
        for (String id : idChildrenHashMap.keySet()) {
            List<String> childrenIds = idChildrenHashMap.get(id);
            Set<String> genesSet = new HashSet<>();
            for (String child_id : childrenIds) {
                List<String> tmp_genes = termToGeneList.get(child_id);
                if (tmp_genes != null) {
                    genesSet.addAll(tmp_genes);
                }
            }
            childrenGenesMap.put(id, genesSet);
        }

        // 基因名称列表
        geneList = new ArrayList<>(Files.readAllLines(Paths.get("temp/" + gene_type + "/geneList.txt")));
        int size = geneList.size();
        // 基因名称到索引的映射
        if (gene_type.equals("Human")) {
            // 根据映射文件获取标识符与基因的映射关系
            HashMap<String, String> idGeneHashMap = new HashMap<>();
            for (String line : Files.readAllLines(Paths.get("buf/Human/id_symbol.map"))) {
                String[] values = line.split("\t");
                idGeneHashMap.put(values[0], values[1]);
            }
            // 初始化基因的标识符到索引的映射
            geneIndexMap = new HashMap<>();
            for (int i = 0; i < size; i++) {
                geneIndexMap.put(idGeneHashMap.get(geneList.get(i)), i);
            }
        } else {
            // 初始化基因名称（标识符）到索引的映射
            geneIndexMap = new HashMap<>();
            for (int i = 0; i < size; i++) {
                geneIndexMap.put(geneList.get(i), i);
            }
        }
        // 初始化基因相关性得分矩阵
        geneSimilarityMatrix = new double[size][size];
        // 读取矩阵数据并填充二维数组
        List<String> lines = Files.readAllLines(Paths.get("result/" + gene_type + "/result" + alpha + ".mat"));
        for (int i = 0; i < lines.size(); i++) {
            String[] values = lines.get(i).split("\t");
            for (int j = 0; j < values.length; j++) {
                double value = Double.parseDouble(values[j]);
                geneSimilarityMatrix[i][j] = value; // 填充矩阵
            }
        }

        // 确保相关性矩阵数据初始化
        isInit = true;
    }

    // 基因名称到索引的映射
    // 数组去重
    static void arrayDistinct(ArrayList<String> arr) {
        HashSet<String> set = new HashSet<>(arr);
        arr.clear();
        arr.addAll(set);
    }

    // 计算基因集距离D
    private static double get_D(String term1, String term2, String gene1, String gene2) {
        // 获取注释基因集合并转换为Set以去重
        HashSet<String> set1 = new HashSet<>(termToGeneList.getOrDefault(term1, new ArrayList<>()));
        HashSet<String> set2 = new HashSet<>(termToGeneList.getOrDefault(term2, new ArrayList<>()));
        // 移除需要比较的基因
        set1.remove(gene1);
        set1.remove(gene2);
        set2.remove(gene1);
        set2.remove(gene2);
        // 计算距离
        double dij_12 = get_dij(set1, set2);
        double dij_21 = get_dij(set2, set1);
        // 计算并集
        set1.addAll(set2);
        return (dij_12 + dij_21) / (2 * set1.size() - dij_12 - dij_21);
    }

    // 计算基因间的距离dij
    private static double get_dij(HashSet<String> set1, HashSet<String> set2) {
        double sum = 0.0;  // 用于存储距离的总和
        double similarity;
        for (String gene1 : set1) {
            double value = 1.0;  // 用于计算每个gene1与set2中所有基因的乘积
            int index1 = geneIndexMap.getOrDefault(gene1, -1); // 获取gene1的索引
            for (String gene2 : set2) {
                if (gene1.equals(gene2)) {
                    value = 0.0; // 两个基因相同，距离为0，接下来无论如何遍历value都是0
                    break;
                }
                int index2 = geneIndexMap.getOrDefault(gene2, -1); // 获取gene2的索引
                if (index1 == -1 || index2 == -1) { // 存在不存在的基因
                    continue;
                }
                // 获取gene1和gene2间的距离：dij = 1 - sim(gene1, gene2)
                similarity = geneSimilarityMatrix[index1][index2];
                value *= (1.0D - similarity);
            }
            sum += value; // 添加当前gene1的乘积到总和
        }
        return sum;
    }

    // 计算路径约束注释信息U
    private static int get_U(String term1, String term2, String parent, String gene1, String gene2) {
        Set<String> genes = new HashSet<>();

        // 遍历术语直接添加基因，避免多个集合的使用
        addGenesFromTerms(new String[]{term1, term2, parent}, genes);

        // 添加术语节点到祖先节点的路径对应基因
        addGenesFromPath(term1, parent, genes);
        addGenesFromPath(term2, parent, genes);

        // 移去比较的两个基因
        genes.remove(gene1);
        genes.remove(gene2);

        return genes.size();
    }

    // 从术语直接添加基因
    private static void addGenesFromTerms(String[] terms, Set<String> genes) {
        for (String term : terms) {
            ArrayList<String> geneList = termToGeneList.get(term);
            if (geneList != null) {
                genes.addAll(geneList);
            }
        }
    }

    // 从路径中获取基因
    private static void addGenesFromPath(String term, String parent, Set<String> genes) {
        // 遍历路径节点，避免使用临时集合
        for (String node : getPathWayTermNode(term, parent)) {
            ArrayList<String> geneList = termToGeneList.get(node);
            if (geneList != null) {
                genes.addAll(geneList);
            }
        }
    }

    // 迭代获取路径上的基因本体术语节点（深度优先搜索）
    private static List<String> getPathWayTermNode(String child, String parent) {
        List<String> nodes = new ArrayList<>();
        if (child.equals(parent)) {
            nodes.add(child);
            return nodes;
        }

        Set<String> visited = new HashSet<>();
        Stack<String> stack = new Stack<>();
        stack.push(parent);
        while (!stack.isEmpty()) {
            String current = stack.pop();

            if (visited.add(current)) { // 标记为已访问
                // 获取当前节点的子节点列表
                List<String> children = idSonHashMap.getOrDefault(current, new ArrayList<>());
                // 遍历子节点
                for (String node : children) {
                    // 获取当前节点对应的子集合，使用 HashSet 提高查找效率
                    Set<String> childSet = new HashSet<>(idChildrenHashMap.getOrDefault(node, new ArrayList<>()));
                    // 检查是否有与当前节点相连的子节点
                    if (childSet.contains(child)) {
                        nodes.add(node);
                        stack.push(node); // 将当前节点加入栈中进行后续处理
                    }
                }
            }
        }
        return nodes;
    }

    // 计算两个术语间的相似性（RWRSM）
    static double get_similarity_of_terms(String term1, String term2, String gene1, String gene2) {
        if (!isSameNamespace(term1, term2)) {
            return 0.0;
        }
        if (term1.equals(term2)) {
            return 1.0;
        }

        // 缓存名称空间
        List<Annotation> annotations = idAnnotationHashMap.get(term1);
        String nameSpace = (annotations != null && !annotations.isEmpty()) ? annotations.get(0).nameSpace : "";
        int count = 0;
        for (Annotation annotation : annotationList) {
            if (annotation.nameSpace.equals(nameSpace)) {
                count++;
            }
        }

        // 计算 D(t1, t2)
        double dab = get_D(term1, term2, gene1, gene2);
        double dab_2 = dab * dab;

        Set<String> ga = new HashSet<>(termToGeneList.get(term1));
        Set<String> gb = new HashSet<>(termToGeneList.get(term2));
        int gaSize = ga.size();
        int gbSize = gb.size();
        // 计算 h(t1, t2)
        double h = (dab_2 * count) + ((1 - dab_2) * Math.max(gaSize, gbSize));

        Set<String> commonParents = new HashSet<>(idParentsHashMap.get(term1));
        commonParents.retainAll(idParentsHashMap.get(term2));
        double maxSimilarity = 0.0;
        for (String id : commonParents) {
            int u = get_U(term1, term2, id, gene1, gene2);
            double f = dab_2 * u + (1 - dab_2) * Math.sqrt(gaSize * gbSize);
            // 直接从预处理结果中获取基因
            Set<String> p_genes = childrenGenesMap.get(id);
            double gp = (p_genes != null) ? p_genes.size() : 0;
            // 计算相似度
            double similarity = calculateSimilarity(count, f, h, gaSize, gbSize, gp);
            maxSimilarity = Math.max(maxSimilarity, similarity);
        }
        return maxSimilarity;
    }

    // 命名空间是否相同
    private static boolean isSameNamespace(String term1, String term2) {
        return idTermHashMap.get(term1) != null && idTermHashMap.get(term2) != null &&
                idTermHashMap.get(term1).namespace.equals(idTermHashMap.get(term2).namespace);
    }

    // 计算术语和术语集合中术语之间的相似度
    private static double calculateSimilarity(int count, double f, double h, int gaSize, int gbSize, double gp) {
        return ((2 * Math.log(count) - 2 * Math.log(f)) / (2 * Math.log(count) - Math.log(gaSize) - Math.log(gbSize))) * (1 - ((h / count) * (gp / count)));
    }

    // Part3
    // 计算术语集合和术语集合中术语之间相似度的最大值
    private static double get_max_similarity_between_goSet_and_goSet(List<String> mainTerms, List<String> compareTerms, HashMap<String, Double> computedPairs, String gene1, String gene2) {
        double maxSum = 0.0;
        for (String t1 : mainTerms) {
            double maxSimilarity = 0.0;
            for (String t2 : compareTerms) {
                String pairKey = t1.compareTo(t2) < 0 ? t1 + "-" + t2 : t2 + "-" + t1; // 使用字符串作为缓存键
                double similarity;
                // 检查缓存
                if (!computedPairs.containsKey(pairKey)) {
                    similarity = Calculator.get_similarity_of_terms(t1, t2, gene1, gene2);
                    computedPairs.put(pairKey, similarity);
                } else {
                    similarity = computedPairs.get(pairKey);
                }
                maxSimilarity = Math.max(maxSimilarity, similarity);
            }
            maxSum += maxSimilarity; // 计算最大相似度求和
        }
        return maxSum;
    }

    // 计算两个基因间的相似性GS**********************************************
    public static double get_similarity_of_genes(String gene1, String gene2) {
        // 得到每个基因对应的术语集合
        ArrayList<String> terms1 = geneToTermList.get(gene1);
        ArrayList<String> terms2 = geneToTermList.get(gene2);
        if (terms1 == null || terms2 == null) {
            return 0;
        }
        // 缓存用于存储计算过的相似度
        double sumSimilarities = 0.0;
        HashMap<String, Double> computedPairs = new HashMap<>();
        // 计算术语和术语集合之间的相似度
        sumSimilarities += get_max_similarity_between_goSet_and_goSet(terms1, terms2, computedPairs, gene1, gene2);
        sumSimilarities += get_max_similarity_between_goSet_and_goSet(terms2, terms1, computedPairs, gene1, gene2);
        // 基因相似度计算
        return sumSimilarities / (terms1.size() + terms2.size());
    }

    // Part4
    // LFC得分计算器
    static class Evaluator {
        // 根据ec编号计算LFC得分的接口*************************************
        public static void LFC(String gene_type, double alpha) throws Exception {
            final Path path = Paths.get("result/" + gene_type + "/lfc" + alpha + ".txt");
            createFileIfNotExists(path);
            init(gene_type, alpha);
            // 使用 TRUNCATE_EXISTING 确保覆盖原有内容
            Files.write(path, new byte[0], StandardOpenOption.TRUNCATE_EXISTING); // 先清空文件内容
            // 计算LFC得分
            for (String ec : ecNumGeneHashMap.keySet()) {
                double lfc = get_lfc(ec);
                Files.write(path, (ec + "\t" + lfc + "\n").getBytes(),
                        StandardOpenOption.APPEND);
            }
        }

        // 相关初始化
        private static void init(String gene_type, double alpha) throws Exception {
            // 基因证据ec编号 哈希表 初始化
            switch (gene_type) {
                case "Yeast":
                case "Arabidopsis_thaliana":
                    // 基因相似度结果 哈希表 初始化
                    similarityResult = Reader.getSimilarityResult(gene_type, alpha);
                    break;
                case "Human":
                    // 基因相似度结果 哈希表 初始化
                    similarityResult = Reader.getHumanSimilarityResult(gene_type, alpha);
                    break;
            }
        }

        // 计算LFC得分
        private static double get_lfc(String ecNumber) {
            // 获取ei集：用于计算LFC组内差异
            HashSet<String> eiSet = ecNumGeneHashMap.get(ecNumber);
            // 获取不相交ej集：用于计算LFC组间差异
            ArrayList<HashSet<String>> ecSet = getNotInterSet(eiSet);
            double value = 0.0D;
            for (HashSet<String> ejSet : ecSet) {
                double sum = 0.0D;
                for (String gene : eiSet) {
                    sum += get_diff(gene, eiSet, ejSet);
                }
                value += sum / eiSet.size();
            }
            double num_ec = ecNumGeneHashMap.size(); // |EC|
            return value / Math.max(num_ec, 1); // 防止除以0
        }

        // 获取与ei不相交的基因集合ej
        private static ArrayList<HashSet<String>> getNotInterSet(HashSet<String> eiSet) {
            ArrayList<HashSet<String>> result = new ArrayList<>();
            for (String ec : ecNumGeneHashMap.keySet()) {
                HashSet<String> ejSet = ecNumGeneHashMap.get(ec);
                boolean flag = true;
                for (String gene : ejSet) {
                    if (eiSet.contains(gene)) {
                        flag = false;
                        break;
                    }
                }
                if (flag) {
                    result.add(ejSet);
                }
            }
            return result;
        }

        // 计算差异对数diff
        private static double get_diff(String gene, HashSet<String> ei, HashSet<String> ej) {
            // 拉普拉斯平滑参数
            double c = 1.0E-10D;
            HashSet<String> eiCopy = new HashSet<>(ei); // 创建ei的副本
            eiCopy.remove(gene);
            //分子
            double top_value = calculateValue(gene, ej, eiCopy.size(), c);
            // 分母
            double bottom_value = calculateValue(gene, eiCopy, ej.size(), c);

            // 确保 bottom_value大于0，以防止 log 计算出错
            if (bottom_value == 0.0) {
                return 0.0;
            }
            return Math.log(top_value / bottom_value);
        }

        // 提取相似度计算的方法
        private static double calculateValue(String gene, HashSet<String> set, int sizeFactor, double c) {
            double value = 0.0;
            for (String g : set) {
                value += (1 - getSimilarity(gene, g) + c);
            }
            return value * sizeFactor;
        }

        // 获取相似度，考虑顺序对称性
        private static double getSimilarity(String gene, String otherGene) {
            Double similarityValue = similarityResult.get(gene + ":" + otherGene);
            if (similarityValue == null) {
                similarityValue = similarityResult.get(otherGene + ":" + gene);
            }
            if (similarityValue.isNaN()) {
                return 0.0;
            }
            return Math.max(similarityValue, 0.0); // 默认返回0.0
        }

    }

}
