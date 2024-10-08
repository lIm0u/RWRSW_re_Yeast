import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Stream;

public class Reader {
    // 读取基因本体文件，获得基因本体数据数组
    public static ArrayList<Term> readOboFile(String obo) throws IOException {
        ArrayList<Term> termList = new ArrayList<>();
        String content = new String(Files.readAllBytes(Paths.get(obo)));
        Pattern termPattern = Pattern.compile("\\[Term][\\s\\S.+?]+?\\n\\n");
        Pattern idPattern = Pattern.compile("(?<=id: )GO:\\d{7}");
        Pattern namePattern = Pattern.compile("(?<=name: ).+?\\\\n");
        Pattern namespacePattern = Pattern.compile("(?<=namespace: ).+?\\n");
        Pattern isaPattern = Pattern.compile("(?<=is_a: )GO:\\d{7}");
        Pattern partofPattern = Pattern.compile("(?<=relationship: part_of )GO:\\d{7}");
        Pattern obsoletePattern = Pattern.compile("is_obsolete: true");
        Matcher m = termPattern.matcher(content);
        while (m.find()) {
            Term termNode = new Term();
            String term = m.group();

            Matcher idm = idPattern.matcher(term);
            if (idm.find()) {
                termNode.id = idm.group().replace(":", "");
            }
            Matcher names = namePattern.matcher(term);
            if (names.find()) {
                termNode.name = names.group().replace(":", "");
            }
            Matcher nsm = namespacePattern.matcher(term);
            if (nsm.find()) {
                termNode.namespace = nsm.group().replace(":", "");
            }
            Matcher isam = isaPattern.matcher(term);
            while (isam.find()) {
                termNode.isID.add(isam.group().replace(":", ""));
            }
            Matcher part_ofm = partofPattern.matcher(term);
            while (part_ofm.find()) {
                termNode.partID.add(part_ofm.group().replace(":", ""));
            }
            termNode.isObsolete = obsoletePattern.matcher(term).find();
            termList.add(termNode);
        }

        return termList;
    }

    // 读取注释信息文件，获得注释信息数据数组
    public static ArrayList<Annotation> readGafFile(String gaf) throws IOException {
        ArrayList<Annotation> annotationList = new ArrayList<>();

        try (BufferedReader br = new BufferedReader(new FileReader(gaf))) {
            String line;
            while ((line = br.readLine()) != null) {
                // 跳过以!开头的注释行
                if (line.startsWith("!")) {
                    continue;
                }
                String[] row = line.split("\t");
                if (row.length >= 12) {
                    Annotation annotationNode = new Annotation();
                    // 第1列：标识符ID
                    annotationNode.geneID = row[2];
                    // 第2列：基因本体术语ID
                    annotationNode.goID = row[4].replace(":", "");
                    // 第3列：证据编码
                    annotationNode.evidenceCode = row[6];
                    // 第4列：基因本体术语所属分支
                    annotationNode.nameSpace = row[8];
                    // 第5列：基因本体术语的同义词
                    String synonymData = row[9];
                    annotationNode.synonym = new ArrayList<>(); // 初始化同义词列表
                    if (synonymData != null && !synonymData.isEmpty()) { // 检查是否不是空
                        String[] synonyms = synonymData.split("\\|");
                        ArrayList<String> synonymList = new ArrayList<>();
                        for (String synonym : synonyms) {
                            // 过滤无效同义词（如不为空或不只包含空白字符）
                            if (synonym != null && !synonym.trim().isEmpty()) {
                                synonymList.add(synonym.trim());
                            }
                        }
                        annotationNode.synonym = synonymList;
                    }
                    annotationList.add(annotationNode);
                }
            }
        }
        return annotationList;
    }

    // 读取基因功能网络文件，获得基因网络节点数据数组
    public static ArrayList<FunctionNet> readNetFile(String net) throws IOException {
        ArrayList<FunctionNet> netList = new ArrayList<>();
        for (String line : Files.readAllLines(Paths.get(net))) {
            FunctionNet netNode = new FunctionNet();
            netNode.gene1 = line.split("\t")[0];
            netNode.gene2 = line.split("\t")[1];
            netNode.value = Double.parseDouble(line.split("\t")[2]);

            netList.add(netNode);
        }
        return netList;
    }

    // 读取基因相似度文件，获得基因相似度数据数组 哈希表
    public static HashMap<String, Double> getSimilarityResult(String gene_type, double alpha) throws IOException {
        HashMap<String, Double> similarityResult = new HashMap<>();
        for (String line : Files.readAllLines(Paths.get("result/" + gene_type + "/similarityResult" + alpha + ".txt"))) {
            String gene1 = line.split("\t")[0];
            String gene2 = line.split("\t")[1];
            Double value = Double.parseDouble(line.split("\t")[2]);
            similarityResult.put(gene1 + ":" + gene2, value);
        }
        return similarityResult;
    }

    // 读取Yeast ec文件，获得ec数据
    public static HashMap<String, HashSet<String>> getYeastEcNumGeneHashMap(ArrayList<String> filePaths) throws IOException {
        // 基因证据ec编号 哈希表 初始化
        HashMap<String, HashSet<String>> ecNumGeneHashMap = new HashMap<>();
        try (Stream<String> lines = Files.lines(Paths.get(filePaths.get(3)))) {
            lines.map(line -> line.split("\t")) // 使用制表符作为分隔符读取文件
                    .filter(values -> values.length >= 4) // 确保有足够的列
                    .forEach(values -> {
                        String ec = values[2];
                        String gene = values[3];
                        ecNumGeneHashMap.computeIfAbsent(ec, k -> new HashSet<>()).add(gene); // 自动处理缺失项
                    });
        } catch (IOException e) {
            throw new IOException(e);
        }
        // 过滤掉只有一个基因的ec编号
        ecNumGeneHashMap.entrySet().removeIf(entry -> entry.getValue().size() <= 1);
        // 移除空字符串的键（如果存在）
        ecNumGeneHashMap.remove("");
        return ecNumGeneHashMap;
    }

    // 读取Arabidopsis thaliana ec文件，获得ec数据
    public static HashMap<String, HashSet<String>> getArabidopsisEcNumGeneHashMap(ArrayList<String> filePaths) throws IOException {
        // 基因证据ec编号 哈希表 初始化
        HashMap<String, HashSet<String>> ecNumGeneHashMap = new HashMap<>();
        try (BufferedReader br = new BufferedReader(new FileReader(filePaths.get(3)))) {
            String line;
            // 如果有则跳过表头
            br.readLine();
            while ((line = br.readLine()) != null) {
                String[] columns = line.split("\t"); // 使用制表符作为分隔符读取文件
                if (columns.length < 8) {
                    continue;
                } // 确保行中有足够的列
                String ecNumber = columns[3].replace("EC-", "").trim(); // EC编号
                String geneName = columns[7].trim(); // 基因名字
                // 检查EC编号和基因名字是否为空
                if (!ecNumber.isEmpty() && !geneName.isEmpty()) {
                    ecNumGeneHashMap.putIfAbsent(ecNumber, new HashSet<>()); // 如果没有则创建新的集合
                    ecNumGeneHashMap.get(ecNumber).add(geneName); // 添加基因名字到相应的EC编号
                }
            }
        } catch (IOException e) {
            throw new IOException(e);
        }
        // 过滤掉只包含一个基因的EC编号
        ecNumGeneHashMap.values().removeIf(entry -> entry.size() <= 1);
        // 移除空字符串的键（如果存在）
        ecNumGeneHashMap.remove("");
        return ecNumGeneHashMap;
    }

    // 读取Human ec文件，获得ec数据
    public static HashMap<String, HashSet<String>> getHumanEcNumGeneHashMap(ArrayList<String> filePaths) throws IOException {
        HashMap<String, HashSet<String>> ecNumGeneHashMap = new HashMap<>();
        try (BufferedReader br = new BufferedReader(new FileReader(filePaths.get(3)))) {
            String line;
            // 跳过表头（如果有）
            br.readLine();
            while ((line = br.readLine()) != null) {
                String[] columns = line.split(","); // 使用逗号作为分隔符
                if (columns.length < 2) {
                    continue; // 确保行中有足够的列
                }
                String ecStr = columns[0].trim(); // 获取ecs列（假设在第一列）
                String geneIdStr = columns[1].trim(); // 获取gene_ids列（假设在第二列）
                if (!ecStr.isEmpty() && !geneIdStr.isEmpty()) {
                    String[] ecList = ecStr.split("\\|"); // 分割ecs
                    String[] geneIdList = geneIdStr.split("\\|"); // 分割gene_ids
                    // 创建对应关系
                    for (int i = 0; i < Math.min(ecList.length, geneIdList.length); i++) {
                        String geneId = geneIdList[i].trim();
                        String ec = ecList[i].trim();
                        ecNumGeneHashMap.computeIfAbsent(geneId, k -> new HashSet<>()).add(ec); // 添加基因ID到相应的EC编号下
                    }
                }
            }
        } catch (IOException e) {
            throw new IOException(e); // 处理IO异常
        }
        // 过滤掉只包含一个基因的EC编号
        ecNumGeneHashMap.entrySet().removeIf(entry -> entry.getValue().size() <= 1);
        // 移除空字符串的键（如果存在）
        ecNumGeneHashMap.remove("");
        return ecNumGeneHashMap;
    }
}
