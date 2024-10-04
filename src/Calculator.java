import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;

import static java.lang.Math.pow;

public class Calculator {
    private static final HashSet<String> buff = new HashSet<>(); // 缓冲区
    static HashMap<String, Double> netHashMap; // 原始网络数据 对应net.txt文件 存储为gene1:gene2 value
    static ArrayList<String> geneList; // 基因列表 本应来自net.txt文件，此处由于进一步筛选因此直接从geneList.txt读入
    static ArrayList<ArrayList<Double>> netMatrix; // 原始网络矩阵,数值来自netHashMap，行列坐标来自geneList
    static HashMap<String, Double> geneSimilarityHashMap; // 随机游走的基因相关性得分矩阵数据哈希表
    static HashMap<String, ArrayList<String>> sonsBuff; // 直接后代子节点缓冲区
    static ArrayList<Annotation> annotationList; // gene.gaf 注释文件数据
    static HashMap<String, ArrayList<String>> termToGeneList; // 本体术语-->术语注释的基因集间的映射表
    static HashMap<String, ArrayList<String>> idParentsHashMap; // 术语网络拓扑图父节点 哈希表: key=本体id value=本体结构信息
    static HashMap<String, ArrayList<String>> idChildrenHashMap; // 术语网络拓扑图后代节点 哈希表: key=本体id value=本体结构信息
    static HashMap<String, ArrayList<String>> idSonHashMap; // 术语网络拓扑图直接儿子节点 哈希表: key=本体id value=本体结构信息
    static HashMap<String, Term> idTermHashMap; // 本体术语项 哈希表: key=基因本体id value=基因本体结构信息
    static HashMap<String, ArrayList<Annotation>> idAnnotationHashMap; // 注释信息数据 哈希表: key=本体ID value=注释信息
    static HashMap<String, ArrayList<String>> geneToTermList;// 基因-->注释基因的本体术语集合间的映射表
    static HashMap<String, HashSet<String>> ecNumGeneHashMap; // 基因证据编号 哈希表: key=证据编码 value=该ec编码的基因集
    static HashMap<String, Double> similarityResult; // 相似度计算结果
    private static boolean isInit = false; // 是否初始化标志

    // Part1
    // 判断文件是否存在，并创建文件
    public static void createFileIfNotExists(Path path) throws IOException {
        if (!Files.exists(path)) {
            Files.createFile(path);
        }
    }

    // 根据参数选择文件读写路径
    static ArrayList<String> getFilePaths(String gene_type) {
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
        createFileIfNotExists(Paths.get("temp/" + gene_type + "/geneList.txt"));
        Files.write(Paths.get("temp/" + gene_type + "/geneList.txt"), sg.toString().getBytes());

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
        createFileIfNotExists(Paths.get("temp/" + gene_type + "/netMatrix.txt"));
        Files.write(Paths.get("temp/" + gene_type + "/netMatrix.txt"), sb.toString().getBytes());
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
        createFileIfNotExists(Paths.get("temp/" + gene_type + "/son.txt"));
        try (BufferedWriter writer = Files.newBufferedWriter(Paths.get("temp/" + gene_type + "/son.txt"))) {
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
        createFileIfNotExists(Paths.get("temp/" + gene_type + "/child.txt"));
        try (BufferedWriter writer = Files.newBufferedWriter(Paths.get("temp/" + gene_type + "/child.txt"))) {
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

    // 相关数据初始化
    public static void data_initializer(String gene_type, double alpha) throws Exception {
        if (isInit) {
            return;
        }
        // 根据参数修改文件路径
        ArrayList<String> filePaths = getFilePaths(gene_type);
        ArrayList<Term> obo_terms = Reader.readOboFile("buf/go-basic.obo");
        // 构建基因本体树相关数据
//        createGOTree(obo_terms, gene_type);

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
            ArrayList<String> childs = new ArrayList<>();
            if (line.split(":").length > 1) {
                childs.addAll(Arrays.asList(line.split(":")[1].split("\t")));
            }
            idChildrenHashMap.put(line.split(":")[0], childs);
        }

        // 术语网络拓扑图直接儿子节点 哈希表 初始化
        idSonHashMap = new HashMap<>();
        for (String line : Files.readAllLines(Paths.get("temp/" + gene_type + "/son.txt"))) {
            ArrayList<String> sons = new ArrayList<>();
            if (line.split(":").length > 1) {
                sons.addAll(Arrays.asList(line.split(":")[1].split("\t")));
            }
            idSonHashMap.put(line.split(":")[0], sons);
        }

        // 注释文件数据 gene.gaf
        annotationList = Reader.readGafFile(filePaths.get(2));

        // 注释信息数据 哈希表 初始化
        idAnnotationHashMap = new HashMap<>();
        for (Annotation annotation : annotationList) {
            idAnnotationHashMap.putIfAbsent(annotation.goID, new ArrayList<>());
            idAnnotationHashMap.get(annotation.goID).add(annotation);
        }

        // 基因本体术语与术语注释的基因集间的映射表 初始化
        termToGeneList = new HashMap<>();
        for (Annotation annotation : annotationList) {
            termToGeneList.putIfAbsent(annotation.goID, new ArrayList<>());
            termToGeneList.get(annotation.goID).add(annotation.geneID);
            if (annotation.synonym != null) {
                termToGeneList.get(annotation.goID).addAll(annotation.synonym);
            }
        }
        for (String term : termToGeneList.keySet()) {
            arrayDistinct(termToGeneList.get(term));
        }

        // 基因-基因本体/同义词-基因本体 哈希表 初始化
        geneToTermList = new HashMap<>();
        for (Annotation annotation : annotationList) {
            geneToTermList.putIfAbsent(annotation.geneID, new ArrayList<>());
            geneToTermList.get(annotation.geneID).add(annotation.goID);
            if (annotation.synonym != null) {
                for (String gene : annotation.synonym) {
                    geneToTermList.putIfAbsent(gene, new ArrayList<>());
                    geneToTermList.get(gene).add(annotation.goID);
                }
            }
        }
        for (String gene : geneToTermList.keySet()) {
            arrayDistinct(geneToTermList.get(gene));
        }

        // 基因名称列表
        ArrayList<String> geneNames = new ArrayList<>(Files.readAllLines(Paths.get("temp/" + gene_type + "/geneList.txt")));
        // 基因相关性得分矩阵数据哈希表
        geneSimilarityHashMap = new HashMap<>();
        int i = 0;
        for (String line : Files.readAllLines(Paths.get("result/" + gene_type + "/result" + alpha + ".mat"))) {
            String[] values = line.split("\t");
            for (int j = 0; j < values.length; j++) {
                double value = Double.parseDouble(values[j]);
                if (value != 0.0) {
                    String key = geneNames.get(i) + ":" + geneNames.get(j);
                    geneSimilarityHashMap.put(key, value);
                }
            }
            i++;
        }
        isInit = true;
    }

    // 数组去重
    static void arrayDistinct(ArrayList<String> arr) {
        HashSet<String> set = new HashSet<>(arr);
        arr.clear();
        arr.addAll(set);
    }

    // 计算基因间的距离dij
    private static double get_dij(Set<String> set1, Set<String> set2) {
        double sum = 0.0;  // 用于存储距离的总和
        for (String gene1 : set1) {
            double product = 1.0;  // 用于计算每个gene1与set2中所有基因的乘积
            for (String gene2 : set2) {
                if (gene1.equals(gene2)) {
                    product = 0; // 不同基因不能相同，如果相同则距离为0
                    continue; // 下一个基因
                }
                // 计算gene1和gene2间的距离
                double similarity = geneSimilarityHashMap.getOrDefault(gene1 + ":" + gene2, 0.0D);
                product *= (1.0D - similarity); // 乘积
            }
            sum += product; // 添加当前gene1的乘积到总和
        }
        return sum;
    }

    // 计算基因集距离D
    private static double get_D(String term1, String term2, String gene1, String gene2) {
        // 获取注释基因集合并转换为Set以去重
        Set<String> set1 = new HashSet<>(termToGeneList.getOrDefault(term1, new ArrayList<>()));
        Set<String> set2 = new HashSet<>(termToGeneList.getOrDefault(term2, new ArrayList<>()));
        // 移除需要比较的基因
        set1.remove(gene1);
        set1.remove(gene2);
        set2.remove(gene1);
        set2.remove(gene2);
        // 计算距离
        double dij_12 = get_dij(set1, set2);
        double dij_21 = get_dij(set2, set1);
        // 计算并集的大小
        set1.addAll(set2);
        // 计算D
        return (dij_12 + dij_21) / (2 * set1.size() - dij_12 - dij_21);
    }

    // 获取路径上的基因本体术语节点
    private static ArrayList<String> getPathWayTermNode(String child, String parent) {
        ArrayList<String> nodes = new ArrayList<>();
        if (child.equals(parent)) {
            nodes.add(child);
            return nodes;
        }
        Set<String> visited = new HashSet<>(); // 使用集合避免重复访问
        getPathWayTermNodeRecursive(child, parent, nodes, visited);
        return nodes;
    }

    // 递归获取路径上的基因本体术语节点
    private static void getPathWayTermNodeRecursive(String child, String parent, List<String> nodes, Set<String> visited) {
        if (!visited.add(parent)) {
            return; // 使用 add() 方法，如果已经存在则返回 false
        }// 标记为已访问
        for (String node : idSonHashMap.getOrDefault(parent, new ArrayList<>())) {
            List<String> children = idChildrenHashMap.get(node);
            if (children != null && children.contains(child)) {
                nodes.add(node);
                getPathWayTermNodeRecursive(child, node, nodes, visited);
            }
        }
    }

    // 从路径中获取基因
    private static void getGenesFromPath(String term, String parent, Set<String> genes) {
        for (String node : getPathWayTermNode(term, parent)) {
            ArrayList<String> geneList = termToGeneList.get(node); // 直接进行 null 检查
            if (geneList != null) {
                genes.addAll(geneList);
            }
        }
    }

    // 计算路径约束注释信息U
    private static int get_U(String term1, String term2, String parent, String gene1, String gene2) {
        Set<String> genes = new HashSet<>();
        // 使用一个集合直接添加所有基因
        for (String term : new String[]{term1, term2, parent}) {
            ArrayList<String> geneList = termToGeneList.get(term);
            if (geneList != null) {
                genes.addAll(geneList);
            }
        }
        // 添加术语节点a到祖先节点的路径对应基因
        getGenesFromPath(term1, parent, genes);
        // 添加术语节点b到祖先节点的路径对应基因
        getGenesFromPath(term2, parent, genes);
        // 移去比较的两个基因
        genes.remove(gene1);
        genes.remove(gene2);
        return genes.size();
    }

    // 计算两个术语间的相似性（RWRSM）
    static double get_similarity_of_terms(String term1, String term2, String gene1, String gene2) {
        if (!isSameNamespace(term1, term2)) {
            return 0.0;
        }
        if (term1.equals(term2)) {
            return 1.0;
        }

        Set<String> commonParents = new HashSet<>(idParentsHashMap.get(term1));
        commonParents.retainAll(idParentsHashMap.get(term2));

        Set<String> ga = new HashSet<>(termToGeneList.get(term1));
        Set<String> gb = new HashSet<>(termToGeneList.get(term2));

        // 缓存名称空间
        String nameSpace = idAnnotationHashMap.get(term1) != null ? idAnnotationHashMap.get(term1).get(0).nameSpace : "";
        int count = (int) annotationList.stream()
                .filter(annotation -> annotation.nameSpace.equals(nameSpace))
                .count();

        double dab = get_D(term1, term2, gene1, gene2);
        double dab_2 = dab * dab;
        // 计算 h(t1, t2)
        double h = (dab_2 * count) + ((1 - dab_2) * Math.max(ga.size(), gb.size()));
        double maxSimilarity = 0.0;
        for (String id : commonParents) {
            int u = get_U(term1, term2, id, gene1, gene2);
            double f = dab_2 * u + (1 - dab_2) * Math.sqrt(ga.size() * gb.size());
            double gp = 0;
            ArrayList<String> p_genes = new ArrayList<>();
            ArrayList<String> childrenIds = idChildrenHashMap.get(id);
            if (childrenIds != null) {
                for (String child_id : childrenIds) {
                    ArrayList<String> tmp_genes = termToGeneList.get(child_id);
                    if (tmp_genes != null) {
                        p_genes.addAll(tmp_genes);
                    }
                }
                arrayDistinct(p_genes);
                gp = p_genes.size();
            }
            double similarity = calculateSimilarity(count, f, h, ga.size(), gb.size(), gp);
            maxSimilarity = Math.max(maxSimilarity, similarity);
        }
        if (maxSimilarity > 1.0) {
            maxSimilarity = 1.0;
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
    private static double get_max_similarity_between_goSet_and_goSet(List<String> mainTerms, List<String> compareTerms, Map<String, Double> computedPairs, String gene1, String gene2) {
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
        Map<String, Double> computedPairs = new ConcurrentHashMap<>();
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
            ArrayList<String> filePaths = getFilePaths(gene_type);
            final Path path = Paths.get("result/" + gene_type + "/lfc" + alpha + ".txt");
            createFileIfNotExists(path);
            init(gene_type, filePaths, alpha);
            // 使用 TRUNCATE_EXISTING 确保覆盖原有内容
            Files.write(path, new byte[0], StandardOpenOption.TRUNCATE_EXISTING); // 先清空文件内容
            // 计算LFC得分
            for (String ec : ecNumGeneHashMap.keySet()) {
                double lfc = get_lfc(ec);
                Files.write(path, (ec + "\t\t" + lfc + "\n").getBytes(),
                        StandardOpenOption.APPEND);
            }
        }

        // 获取只有ei才有的基因集合ej
        private static ArrayList<HashSet<String>> getNotInterSet(HashSet<String> eiSet) {
            ArrayList<HashSet<String>> result = new ArrayList<>();
            for (String ec : ecNumGeneHashMap.keySet()) {
                HashSet<String> ej = ecNumGeneHashMap.get(ec);
                boolean flag = true;
                for (String gene : ej) {
                    if (eiSet.contains(gene)) {
                        flag = false;
                        break;
                    }
                }
                if (flag) {
                    result.add(ej);
                }
            }
            return result;
        }

        // 计算差异对数diff
        private static double get_diff(String gene, HashSet<String> ei, HashSet<String> ej) {
            double top_value = calculateValue(gene, ej, ei.size());
            double bottom_value = calculateValue(gene, ei, ej.size());

            // 确保 top_value 和 bottom_value 都大于0，以防止 log 计算出错
            top_value = Math.max(top_value, 1.0E-10);
            bottom_value = Math.max(bottom_value, 1.0E-10);
            if (bottom_value == 0.0) {
                return 0.0;
            }
            // if (result < 0) {
//                System.out.print(" < 0: " + gene + " " + ei + " / " + ej + top_value + " / " + bottom_value);
//                return 0;
//            }
//            if (result > 1) {
//                System.out.print(" > 1: " + gene + " " + ei + " / " + ej + "\n");
//                return 1;
//            }
            return Math.log(top_value / bottom_value);
        }

        // 提取相似度计算的方法
        private static double calculateValue(String gene, HashSet<String> set, int sizeFactor) {
            double value = 0.0;
            for (String g : set) {
                value += (1 - getSimilarity(gene, g) + 1.0E-10);
            }
            return value * sizeFactor;
        }

        // 获取相似度，考虑顺序对称性
        private static double getSimilarity(String gene, String otherGene) {
            Double similarityValue = similarityResult.get(gene + ":" + otherGene);
            if (similarityValue == null) {
                similarityValue = similarityResult.get(otherGene + ":" + gene);
            }
            return similarityValue != null ? Math.max(similarityValue, 0.0) : 0.0; // 默认返回0.0
        }

        // 计算LFC得分
        private static double get_lfc(String ecNumber) {
            HashSet<String> ei = ecNumGeneHashMap.get(ecNumber);
            // 获取不相交ej集
            ArrayList<HashSet<String>> ecSet = getNotInterSet(ei);
            double value = 0.0D;
            for (HashSet<String> ej : ecSet) {
                double sum = 0.0D;
                for (String gene : ei) {
                    sum += get_diff(gene, ei, ej);
                }
                value += sum / ei.size();
            }
            double num_ec = ecNumGeneHashMap.size(); // |EC|
            return value / Math.max(num_ec, 1); // 防止除以0
        }

        // 相关初始化
        private static void init(String gene_type, ArrayList<String> filePaths, double alpha) throws Exception {
            // 基因相似度结果 哈希表 初始化
            similarityResult = Reader.getSimilarityResult(gene_type, alpha);
            // 基因证据ec编号 哈希表 初始化
            switch (gene_type) {
                case "Yeast":
                    ecNumGeneHashMap = Reader.getYeastEcNumGeneHashMap(filePaths);
                    break;
                case "Human":
                    ecNumGeneHashMap = Reader.getHumanEcNumGeneHashMap(filePaths);
                    break;
                case "Arabidopsis_thaliana":
                    ecNumGeneHashMap = Reader.getArabidopsisEcNumGeneHashMap(filePaths);
                    break;
            }
        }
    }

}
