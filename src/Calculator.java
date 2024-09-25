import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.OpenOption;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.*;

public class Calculator {
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
    private static HashSet<String> buff = new HashSet<>(); // 缓冲区
    private static boolean isInit = false; // 是否初始化标志

    // Part1
    // 计算并存储原始基因网络数据矩阵 netMatrix
    private static void getNetMatrix() throws IOException {
        HashSet<String> tmp = new HashSet<>();

        // 初始化 netHashMap
        netHashMap = new HashMap<>();
        for (FunctionNet edge : Reader.readNetFile("data/net.txt")) {
            netHashMap.put(edge.gene1 + ":" + edge.gene2, edge.value / 10.0D); //做归一化处理
            String gene1 = edge.gene1;
            String gene2 = edge.gene2;
            tmp.add(gene1);
            tmp.add(gene2);
        }
        geneList = new ArrayList<>(tmp);

        // 存储 geneList.txt
        StringBuilder sbg = new StringBuilder();
        for (String line : geneList) {
            sbg.append(line).append("\n");
        }
        if (!Files.exists(Paths.get("temp/geneList.txt")))
            Files.createFile(Paths.get("temp/geneList.txt"));
        Files.write(Paths.get("temp/geneList.txt"), sbg.toString().getBytes());

        // 计算 netMatrix
        netMatrix = new ArrayList<>();
        netMatrix.ensureCapacity(4172);
        for (int i = 0; i < geneList.size(); i++) {
            netMatrix.add(new ArrayList<>());
            netMatrix.get(i).ensureCapacity(4172);
            for (String s : geneList) {
                String key = geneList.get(i) + ":" + s;
                netMatrix.get(i).add(netHashMap.getOrDefault(key, 0.0D));
            }
        }

        // 存储 netMatrix矩阵
        StringBuilder sb = new StringBuilder();
        for (ArrayList<Double> line : netMatrix) {
            for (Iterator<Double> iterator = line.iterator(); iterator.hasNext(); ) {
                double value = iterator.next();
                sb.append(value).append("\t");
            }
            sb.append("\n");
        }
        if (!Files.exists(Paths.get("temp/netMatrix.txt")))
            Files.createFile(Paths.get("temp/netMatrix.txt"));
        Files.write(Paths.get("temp/netMatrix.txt"), sb.toString().getBytes());
    }

    // 计算权重矩阵
    private static ArrayList<ArrayList<Double>> getWMatrix() throws Exception {
        double sita_2 = parameterVariable();
        ArrayList<ArrayList<Double>> wM = new ArrayList<>();
        wM.ensureCapacity(4172);
        for (int i = 0; i < netMatrix.size(); i++) {
            wM.add(new ArrayList<>());
            wM.get(i).ensureCapacity(4172);
            for (int j = 0; j < netMatrix.size(); j++) {
                wM.get(i).add(Math.exp(-(euclideanDistance(netMatrix.get(i), netMatrix.get(j)) / sita_2)));
            }
        }
        return wM;
    }

    // 计算参数变量的平方
    private static double parameterVariable() throws Exception {
        double sum = 0.0D;
        for (int i = 0; i < netMatrix.size(); i++) {
            for (int j = 0; j < netMatrix.size(); j++) {
                if (i != j)
                    sum += euclideanDistance(netMatrix.get(i), netMatrix.get(j));
            }
        }
        sum /= netMatrix.size();
        sum /= (netMatrix.size() - 1);
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
        ArrayList<ArrayList<Double>> pMatrix = new ArrayList<>();
        ArrayList<ArrayList<Double>> wMatrix = getWMatrix();
        pMatrix.ensureCapacity(4172);
        for (int i = 0; i < netMatrix.size(); i++) {
            pMatrix.add(new ArrayList<>());
            pMatrix.get(i).ensureCapacity(4172);
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
    private static void storagePTMatrix() throws Exception {
        ArrayList<ArrayList<Double>> pm = getPMatrix();
        StringBuilder sb = new StringBuilder();
        for (ArrayList<Double> line : pm) {
            for (Iterator<Double> iterator = line.iterator(); iterator.hasNext(); ) {
                double value = iterator.next();
                sb.append(value).append("\t");
            }
            sb.append("\n");
        }
        if (!Files.exists(Paths.get("temp/pMatrix.txt")))
            Files.createFile(Paths.get("temp/pMatrix.txt"));
        Files.write(Paths.get("temp/pMatrix.txt"), sb.toString().getBytes());
    }

    // 计算随机游走矩阵即基因之间的相关性得分矩阵*****************************************
    public static void getRWRMatrix(double arg) throws Exception {
//        getNetMatrix();
//        storagePTMatrix();
//        RandWalk.walk(arg);
    }

    // Part2
    // 获取后代节点
    private static void getChilds(String goID) {
        if (!sonsBuff.get(goID).isEmpty()) {
            for (String go : sonsBuff.get(goID))
                getChilds(go);
        }
        buff.add(goID);
    }

    // 基因本体树构建
    private static void createGOTree() throws Exception {
        // 生成 son.txt
        sonsBuff = new HashMap<>();
        ArrayList<String> temp;
        for (Term term : Reader.readOboFile("data/onto.obo")) {
            if (term.isObsolote) {
                continue;
            }
            if (!sonsBuff.containsKey(term.id))
                sonsBuff.put(term.id, new ArrayList<>());

            temp = term.isID;
            for (String parent : temp) {
                if (!sonsBuff.containsKey(parent))
                    sonsBuff.put(parent, new ArrayList<>());
                sonsBuff.get(parent).add(term.id);
            }
            temp = term.partID;
            for (String parent : temp) {
                if (!sonsBuff.containsKey(parent))
                    sonsBuff.put(parent, new ArrayList<>());
                sonsBuff.get(parent).add(term.id);
            }
        }
        if (!Files.exists(Paths.get("temp/son.txt")))
            Files.createFile(Paths.get("temp/son.txt"));
        for (String parent : sonsBuff.keySet()) {
            StringBuilder sb = new StringBuilder();
            sb.append(parent).append(":");
            for (String child : sonsBuff.get(parent)) {
                sb.append(child).append("\t");
            }
            sb.append("\n");
            Files.write(Paths.get("temp/son.txt"), sb.toString().getBytes(),
                    new OpenOption[]{StandardOpenOption.APPEND});
        }

        // 生成 child.txt
        HashMap<String, ArrayList<String>> childs = new HashMap<>();
        for (String parent : sonsBuff.keySet()) {
            if (!childs.containsKey(parent))
                childs.put(parent, new ArrayList<>());
            temp = sonsBuff.get(parent);
            for (String child : temp) {
                getChilds(child);
            }
            childs.get(parent).addAll(buff);
            buff.clear();
        }
        if (!Files.exists(Paths.get("temp/child.txt")))
            Files.createFile(Paths.get("temp/child.txt"));
        for (String parent : childs.keySet()) {
            StringBuilder sb = new StringBuilder();
            sb.append(parent).append(":");
            for (String child : childs.get(parent)) {
                sb.append(child).append("\t");
            }
            sb.append("\n");
            Files.write(Paths.get("temp/child.txt"), sb.toString().getBytes(),
                    new OpenOption[]{StandardOpenOption.APPEND});
        }
    }

    // 相关数据初始化
    public static void data_initializer() throws Exception {
        if (isInit) {
            return;
        }

        // 构建基因本体树相关数据
//        createGOTree();

        // 本体术语项 哈希表 初始化
        idTermHashMap = new HashMap<>();
        for (Term term : Reader.readOboFile("data/onto.obo")) {
            if (!term.isObsolote) {
                idTermHashMap.put(term.id, term);
            }
        }


        // 术语网络拓扑图父节点 哈希表 初始化
        idParentsHashMap = new HashMap<>();
        for (String id : idTermHashMap.keySet()) {
            ArrayList<String> parents = new ArrayList<>();
            idParentsHashMap.put(id, parents);
            //已经添加的术语
            HashSet<String> added = new HashSet<>();
            //等待添加的术语
            ArrayList<String> adding = new ArrayList<>();
            adding.addAll(idTermHashMap.get(id).partID);
            adding.addAll(idTermHashMap.get(id).isID);
            while (!adding.isEmpty()) {
                try {
                    if (!added.contains(adding.get(0))) {
                        parents.add(adding.get(0));
                        Term t = idTermHashMap.get(adding.get(0));
                        if (t.partID != null && !t.partID.isEmpty())
                            adding.addAll(t.partID);
                        if (t.isID != null && !t.isID.isEmpty())
                            adding.addAll(t.isID);
                    }
                    adding.remove(0);
                } catch (Exception e) {
                    System.out.println(adding.get(0));
                }
            }
        }

        // 术语网络拓扑图后代节点 哈希表 初始化
        idChildrenHashMap = new HashMap<>();
        for (String line : Files.readAllLines(Paths.get("temp/child.txt"))) {
            ArrayList<String> childs = new ArrayList<>();
            if (line.split(":").length > 1) {
                childs.addAll(Arrays.asList(line.split(":")[1].split("\t")));
            }
            idChildrenHashMap.put(line.split(":")[0], childs);
        }

        // 术语网络拓扑图直接儿子节点 哈希表 初始化
        idSonHashMap = new HashMap<>();
        for (String line : Files.readAllLines(Paths.get("temp/son.txt"))) {
            ArrayList<String> sons = new ArrayList<>();
            if (line.split(":").length > 1) {
                sons.addAll(Arrays.asList(line.split(":")[1].split("\t")));
            }
            idSonHashMap.put(line.split(":")[0], sons);
        }

        // 注释文件数据 gene.gaf
        annotationList = Reader.readGafFile("data/gene.gaf");

        // 注释信息数据 哈希表 初始化
        idAnnotationHashMap = new HashMap<>();
        for (Annotation annotation : annotationList) {
            if (!idAnnotationHashMap.containsKey(annotation.goID))
                idAnnotationHashMap.put(annotation.goID, new ArrayList<>());
            idAnnotationHashMap.get(annotation.goID).add(annotation);
        }

        // 基因本体术语与术语注释的基因集间的映射表 初始化
        termToGeneList = new HashMap<>();
        for (Annotation annotation : annotationList) {
            if (!termToGeneList.containsKey(annotation.goID)) {
                termToGeneList.put(annotation.goID, new ArrayList<>());
            }
            termToGeneList.get(annotation.goID).add(annotation.geneID);
            if (annotation.synonym != null) {
                termToGeneList.get(annotation.goID).addAll(annotation.synonym);
            }
        }
        for (String term : termToGeneList.keySet())
            arrayDistinct(termToGeneList.get(term));

        // 基因-基因本体/同义词-基因本体 哈希表 初始化
        geneToTermList = new HashMap<>();
        for (Annotation annotation : annotationList) {
            if (!geneToTermList.containsKey(annotation.geneID)) {
                geneToTermList.put(annotation.geneID, new ArrayList<>());
            }
            geneToTermList.get(annotation.geneID).add(annotation.goID);
            if (annotation.synonym != null) {
                for (String gene : annotation.synonym) {
                    if (!geneToTermList.containsKey(gene)) {
                        geneToTermList.put(gene, new ArrayList<>());
                    }
                    geneToTermList.get(gene).add(annotation.goID);
                }
            }
        }
        for (String gene : geneToTermList.keySet())
            arrayDistinct(geneToTermList.get(gene));

        // 基因名称列表 geneList
        geneList = new ArrayList<>();
        geneList.addAll(Files.readAllLines(Paths.get("temp/goldList.txt")));

        // 基因相关性得分矩阵数据哈希表 geneSimilarityMatrix
        geneSimilarityHashMap = new HashMap<>();
        int i = 0;
        for (String line : Files.readAllLines(Paths.get("result/result.mat"))) {
            int j = 0;
            for (String value : line.split("\t")) {
                if (!value.equals("0.0")) {
                    String key = geneList.get(i) + ":" + geneList.get(j);
                    geneSimilarityHashMap.put(key, Double.parseDouble(value));
                }
                j++;
            }
            i++;
        }
        isInit = true;
    }

    // 数组去重
    private static void arrayDistinct(ArrayList<String> arr) {
        HashSet<String> set = new HashSet<>(arr);
        arr.clear();
        arr.addAll(set);
    }

    // 计算基因间的距离dij
    private static double get_dij(HashSet<String> set1, HashSet<String> set2) {
        double sum = 0.0D;
        for (String gene1 : set1) {
            double value = 1.0D;
            for (String gene2 : set2) {
                String key = gene1 + ":" + gene2;
                if (gene1.equals(gene2)) {
                    value *= 0.0D;
                } else {
                    value *= (1.0D - geneSimilarityHashMap.getOrDefault(key, 0.0D));
                }
            }
            sum += value;
        }
        return sum;
    }

    // 计算基因集距离D
    private static double get_D(String term1, String term2, String gene1, String gene2, HashMap<String, Double> distanceCache) {
        // 确保使用字典序较小的术语作为键的一部分
        String key1 = term1.compareTo(term2) < 0 ? term1 : term2;
        String key2 = term1.compareTo(term2) < 0 ? term2 : term1;

        String distanceKey = key1 + ":" + key2; // 使用统一的距离键
        if (distanceCache.containsKey(distanceKey)) {
            return distanceCache.get(distanceKey);
        }
        // 使用 HashSet 来存储基因集合
        HashSet<String> set1 = new HashSet<>(termToGeneList.get(term1));
        HashSet<String> set2 = new HashSet<>(termToGeneList.get(term2));
        // 移除基因
        set1.remove(gene1);
        set1.remove(gene2);
        set2.remove(gene1);
        set2.remove(gene2);
        HashSet<String> sum_set = new HashSet<>();
        sum_set.addAll(set1);
        sum_set.addAll(set2);
        // 计算距离
        // 计算距离
        double dij_12 = get_dij(set1, set2);
        double dij_21 = get_dij(set2, set1);
        double value = (dij_12 + dij_21) / ((2 * sum_set.size()) - dij_12 - dij_21);
        distanceCache.put(distanceKey, value);
        return value;
    }

    // 获取路径上的基因本体术语节点
    private static ArrayList<String> getPathWayTermNode(String child, String parent) {
        ArrayList<String> nodes = new ArrayList<>();
        if (child.equals(parent)) {
            nodes.add(child);
            return nodes;
        }
        if (idSonHashMap.get(parent) != null) {
            for (String node : idSonHashMap.get(parent)) {
                ArrayList<String> childs = idChildrenHashMap.get(node);
                if (childs != null && childs.contains(child)) {
                    nodes.add(node);
                    nodes.addAll(getPathWayTermNode(child, node));
                }
            }
        }
        return nodes;
    }

    // 计算路径约束注释信息U
    private static int get_U(String term1, String term2, String parent, String gene1, String gene2) {
        HashSet<String> genes = new HashSet<>();
        ArrayList<String> a_genes = termToGeneList.get(term1);
        if (a_genes != null)
            genes.addAll(a_genes);
        ArrayList<String> b_genes = termToGeneList.get(term2);
        if (b_genes != null)
            genes.addAll(b_genes);
        ArrayList<String> p_genes = termToGeneList.get(parent);
        if (p_genes != null)
            genes.addAll(p_genes);
        // 术语节点a到祖先节点的路径
        ArrayList<String> pathAP = getPathWayTermNode(term1, parent);
        for (String node : pathAP) {
            try {
                ArrayList<String> tmp_genes = termToGeneList.get(node);
                if (tmp_genes != null)
                    genes.addAll(tmp_genes);
            } catch (Exception e) {
                System.out.println(e.getMessage());
            }
        }
        // 术语节点b到祖先节点的路径
        ArrayList<String> pathBP = getPathWayTermNode(term2, parent);
        for (String node : pathBP) {
            try {
                ArrayList<String> tmp_genes = termToGeneList.get(node);
                if (tmp_genes != null)
                    genes.addAll(tmp_genes);
            } catch (Exception e) {
                System.out.println(e.getMessage());
            }
        }
        // 移去比较的两个基因
        genes.remove(gene1);
        genes.remove(gene2);
        return genes.size();
    }

    // 计算两个术语间的相似性（RWRSM）**************************************
    private static double get_similarity_of_terms(String term1, String term2, String gene1, String gene2, Double distance) {
        if (idTermHashMap.get(term1) == null || idTermHashMap.get(term2) == null)
            return 0.0D;
        if (!(idTermHashMap.get(term1)).namespace.equals((idTermHashMap.get(term2)).namespace))
            return 0.0D;
        // 相关参数
        ArrayList<String> p = new ArrayList<>();
        for (String id : idParentsHashMap.get(term1)) {
            if (idParentsHashMap.get(term2).contains(id))
                p.add(id);
        }
        ArrayList<String> ga = termToGeneList.get(term1);
        ArrayList<String> gb = termToGeneList.get(term2);
        int count = (int) annotationList.stream().filter(annotation -> annotation.nameSpace.equals(((Annotation) ((ArrayList) idAnnotationHashMap.get(term1)).get(0)).nameSpace)).count();
        // 计算距离
        double dab = distance; // 改成直接读取
        double dab_2 = dab * dab;
        double maxSimilarity = 0.0D;
        // 计算相似度
        arrayDistinct(p);
        for (String id : p) {
            int u = get_U(term1, term2, id, gene1, gene2);
            double f = dab_2 * u + (1.0D - dab_2) * Math.sqrt((ga.size() * gb.size()));
            int max_gab = Math.max(ga.size(), gb.size());
            double h = dab_2 * count + (1.0D - dab_2) * max_gab;
            int gp = 0;
            ArrayList<String> p_genes = new ArrayList<>();
            if (idChildrenHashMap.get(id) != null) {
                for (String child_id : idChildrenHashMap.get(id)) {
                    ArrayList<String> tmp_genes = termToGeneList.get(child_id);
                    if (tmp_genes != null)
                        p_genes.addAll(tmp_genes);
                }
                arrayDistinct(p_genes);
                gp = p_genes.size();
            }
            double similarity = (2.0D * Math.log(count) - 2.0D * Math.log(f)) / (2.0D * Math.log(count) - Math.log(ga.size()) - Math.log(gb.size())) * (1.0D - h * gp / (count * count));
            maxSimilarity = Math.max(maxSimilarity, similarity);
        }

        return maxSimilarity;
    }

    // Part3
    // 计算术语和术语集合中术语之间相似度的最大值
    private static double get_max_similarity_between_go_and_goSet(String goID, ArrayList<String> goSet, String gene1,
                                                                  String gene2, HashMap<HashSet<String>, Double> computedPairs, HashMap<String, Double> distanceCache) {
        // 计算所有术语对的基因集距离
        int size = goSet.size();
        double[] distance = new double[size];
        double minDistance = Double.MAX_VALUE;
        double maxDistance = Double.MIN_VALUE;


        // 计算距离并找出最小值和最大值
        for (int i = 0; i < size; i++) {
            distance[i] = get_D(goID, goSet.get(i), gene1, gene2, distanceCache);
            minDistance = Math.min(minDistance, distance[i]);
            maxDistance = Math.max(maxDistance, distance[i]);
        }

        // 进行归一化并计算最大相似度
        double maxSimilarity = 0.0D;
        if (maxDistance > minDistance) { // 避免除零错误
            for (int i = 0; i < size; i++) {
                // 归一化
                double normalizedValue = (distance[i] - minDistance) / (maxDistance - minDistance);
                String currentTerm = goSet.get(i);
                HashSet<String> pairKey = new HashSet<>(Arrays.asList(goID, currentTerm));
                double tmpSimilarity;

                // 检查是否计算过
                if (!computedPairs.containsKey(pairKey)) {
                    tmpSimilarity = get_similarity_of_terms(goID, currentTerm, gene1, gene2, normalizedValue);
                    computedPairs.put(pairKey, tmpSimilarity); // 记录计算过的相似度
                } else {
                    tmpSimilarity = computedPairs.get(pairKey); // 使用缓存的结果
                }
                maxSimilarity = Math.max(maxSimilarity, tmpSimilarity);
            }
        } else {
            maxSimilarity = 0.0D;
        }
        return maxSimilarity;
    }

    // 计算两个基因间的相似性GS**********************************************
    public static double get_similarity_of_genes(String gene1, String gene2) {
        // 得到每个基因对应的术语集合
        ArrayList<String> terms1 = geneToTermList.get(gene1);
        ArrayList<String> terms2 = geneToTermList.get(gene2);
        if (terms1 == null || terms2 == null) {
            return 0.0D;
        }

        // 用于存储计算过的相似度
        HashMap<HashSet<String>, Double> computedPairs = new HashMap<>();
        // 计算距离缓存
        HashMap<String, Double> distanceCache = new HashMap<>();
        double termSimilaritySum = 0.0D;

        // 计算 termSimilaritySum
        for (String id : terms1) {
            termSimilaritySum += get_max_similarity_between_go_and_goSet(id, terms2, gene1, gene2, computedPairs, distanceCache);
        }
        for (String id : terms2) {
            termSimilaritySum += get_max_similarity_between_go_and_goSet(id, terms1, gene1, gene2, computedPairs, distanceCache);
        }

        // 基因相似度计算
        return termSimilaritySum / (terms1.size() + terms2.size());
    }

    // Part4
    // LFC得分计算器
    static class Evaluator {
        // 根据ec编号计算LFC得分的接口*************************************
        public static void LFC() throws Exception {
            init();

            if (!Files.exists(Paths.get("result/lfc.txt")))
                Files.createFile(Paths.get("result/lfc.txt"));

            long s = System.currentTimeMillis();
            int number = 1;
            for (String ec : Calculator.ecNumGeneHashMap.keySet()) {
                System.out.println("开始计算第" + number + "个ec分组的得分");
                long start = System.currentTimeMillis();

                double lfc = get_lfc(ec);
                Files.write(Paths.get("result/lfc.txt"), (ec + "\t\t" + lfc + "\n").getBytes(),
                        new OpenOption[]{StandardOpenOption.APPEND});

                long end = System.currentTimeMillis();
                System.out.println(
                        "已算完第" + number++ + "个ec分组的得分," + "耗时:" + (end - start) + "ms,总耗时:" + (end - s) + "ms");
            }
        }

        // 获取ej
        private static ArrayList<HashSet<String>> get_ej_set(HashSet<String> ei) {
            ArrayList<HashSet<String>> result = new ArrayList<>();
            for (String ec : Calculator.ecNumGeneHashMap.keySet()) {
                HashSet<String> ej = Calculator.ecNumGeneHashMap.get(ec);
                boolean flag = true;
                for (String gene : ej) {
                    if (ei.contains(gene)) {
                        flag = false;
                        break;
                    }
                }
                if (flag)
                    result.add(ej);
            }
            return result;
        }

        // 计算差异对数diff
        private static double get_diff(String gene, HashSet<String> ei, HashSet<String> ej) {
            if (gene == null || ei == null || ej == null) {
                throw new IllegalArgumentException("输入参数不能为null");
            }
            Double value = similarityResult.get(gene + ":" + ej);
            if (value == null) {
                // 处理相似度未找到的情况，例如返回0, 或抛出异常
                return 0.0;
            }
            double top_value = 0.0D;
            double bottom_value = 0.0D;
            double c = 1.0E-10D; // 拉普拉斯平滑参数

            for (String gj : ej) {
                top_value += (1 - similarityResult.get(gene + ":" + gj) + c);
            }
            top_value *= ei.size();
            for (String gi : ei) {
                bottom_value += (1 - similarityResult.get(gene + ":" + gi) + c);
            }
            bottom_value *= ej.size();

            return Math.log(top_value / bottom_value);
        }

        // 计算LFC得分
        private static double get_lfc(String ecNumber) {
            HashSet<String> ei = Calculator.ecNumGeneHashMap.get(ecNumber);
            if (ei == null) {
                // 处理未找到的情况，例如返回默认值，抛出异常等
                return 0.0;
            }
            ArrayList<HashSet<String>> ec = get_ej_set(ei);
            double value = 0.0D;
            for (HashSet<String> ej : ec) {
                double sum = 0.0D;
                for (String gene : ei)
                    sum += get_diff(gene, ei, ej);
                value += sum / ei.size();
            }
            double num_ec = Calculator.ecNumGeneHashMap.size();
            return value / num_ec;
        }

        // 相关初始化
        private static void init() throws Exception {
            // 基因相似度结果 哈希表 初始化
            similarityResult = new HashMap<>();
            for (String line : Files.readAllLines(Paths.get("result/similarityResult.txt"))) {
                String gene1 = line.split("\t")[0];
                String gene2 = line.split("\t")[1];
                Double value = Double.parseDouble(line.split("\t")[2]);
                similarityResult.put(gene1 + ":" + gene2, value);
            }

            // 基因证据编号 哈希表 初始化
            ecNumGeneHashMap = new HashMap<>();
            for (String line : Files.readAllLines(Paths.get("data/ec.tab"))) {
                String[] values = line.split("\t");
                if (values.length < 4)
                    continue;
                String ec = values[2];
                String gene = values[3];
                if (!ecNumGeneHashMap.containsKey(ec))
                    ecNumGeneHashMap.put(ec, new HashSet<>());
                ecNumGeneHashMap.get(ec).add(gene);
            }
            HashMap<String, HashSet<String>> tmp = new HashMap<>();
            for (String ec : ecNumGeneHashMap.keySet()) {
                if (ecNumGeneHashMap.get(ec).size() > 1)
                    tmp.put(ec, ecNumGeneHashMap.get(ec));
            }
            ecNumGeneHashMap.clear();
            ecNumGeneHashMap.putAll(tmp);
            ecNumGeneHashMap.remove("");
        }
    }

}
