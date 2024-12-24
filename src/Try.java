import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

public class Try {
    public static void main(String[] args) {
//        String filePath = "result/Yeast/lfc0.9.txt";
        String filePath = "result/Arabidopsis_thaliana/lfc0.9.txt";
        Map<String, Double> lfcMap = readLfcFile(filePath);
        Map<String, Double> filteredLfcMap = filterLfcMapByThreshold(lfcMap, 22);

        List<Double> scores = new ArrayList<>(filteredLfcMap.values());
        Collections.sort(scores);
        printQuantiles(scores);
//        // 基因相似度结果 哈希表 初始化
//        HashSet<String> geneSet = new HashSet<>();
//        HashMap<String, Double> similarityResult = Reader.getSimilarityResult("Arabidopsis_thaliana", 0.9);
//        for (String genes : similarityResult.keySet()) {
//            String[] geneName = genes.split(":");
//            geneSet.add(geneName[0]);
//            geneSet.add(geneName[1]);
//        }
//        System.out.println("共有基因数量：" + geneSet.size());
    }

    // 读取文件并填充 lfcMap
    private static Map<String, Double> readLfcFile(String filePath) {
        Map<String, Double> lfcMap = new HashMap<>();

        try (BufferedReader br = new BufferedReader(new FileReader(filePath))) {
            String line;
            while ((line = br.readLine()) != null) {
                String[] parts = line.trim().split("\\s+");
                if (parts.length == 2 && !parts[1].equals("NaN") && !parts[1].startsWith("-")) {
                    String ecNumber = parts[0];
                    double lfcScore = Double.parseDouble(parts[1]);
                    lfcMap.put(ecNumber, lfcScore);
                }
            }
        } catch (IOException e) {
            throw new RuntimeException("读取文件时发生错误: " + e.getMessage(), e);
        }
        System.out.println(lfcMap.size());
        return lfcMap;
    }

    // 根据阈值过滤 LFC map
    private static Map<String, Double> filterLfcMapByThreshold(Map<String, Double> lfcMap, double threshold) {
        Map<String, Double> filteredMap = new HashMap<>();
        for (Map.Entry<String, Double> entry : lfcMap.entrySet()) {
            if (entry.getValue() <= threshold) {
                filteredMap.put(entry.getKey(), entry.getValue());
            }
        }
        return filteredMap;
    }

    // 打印分位数
    private static void printQuantiles(List<Double> scores) {
        double q1 = calculateQuantile(scores, 0.25);
        double median = calculateQuantile(scores, 0.5);
        double q3 = calculateQuantile(scores, 0.75);

        System.out.println("第一分位数(Q1): " + q1);
        System.out.println("中位数(Q2): " + median);
        System.out.println("第三分位数(Q3): " + q3);
    }

    // 计算分位数
    private static double calculateQuantile(List<Double> sortedScores, double quantile) {
        if (sortedScores.isEmpty()) {
            return Double.NaN; // 返回NaN以处理空列表
        }

        int n = sortedScores.size();
        double index = quantile * (n - 1);
        int lowerIndex = (int) Math.floor(index);
        int upperIndex = (int) Math.ceil(index);

        if (lowerIndex == upperIndex) {
            return sortedScores.get(lowerIndex);
        } else {
            return (sortedScores.get(lowerIndex) + sortedScores.get(upperIndex)) / 2.0;
        }
    }

}