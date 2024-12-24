import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

public class LfcPair {
    public static HashMap<String, HashSet<String>> getLfc(String gene_type) throws IOException {
        HashMap<String, HashSet<String>> ecNumGeneHashMap;

        // 基因证据ec编号 哈希表 初始化
        switch (gene_type) {
            case "Yeast":
                ecNumGeneHashMap = Reader.getYeastEcNumGeneHashMap(Calculator.getFilePaths(gene_type));
                break;
            case "Human":
                ecNumGeneHashMap = Reader.getHumanEcNumGeneHashMap(Calculator.getFilePaths(gene_type));
                break;
            case "Arabidopsis_thaliana":
                ecNumGeneHashMap = Reader.getArabidopsisEcNumGeneHashMap(Calculator.getFilePaths(gene_type));
                break;
            default:
                ecNumGeneHashMap = new HashMap<>();
                break;
        }

        // 创建一个HashSet来存储所有唯一的基因
        HashSet<String> uniqueGenes = new HashSet<>();
        ecNumGeneHashMap.values().forEach(uniqueGenes::addAll);

        // 转换为List以便后续处理
        List<String> uniqueGenesList = new ArrayList<>(uniqueGenes);

        if (gene_type.equals("Human")) {
            // 通过映射文件将基因名转换为标识符
            uniqueGenesList = Reader.getHumanGeneNameToIdentifierMap(uniqueGenesList);
        }

        // 确保文件路径存在
        final Path path = Paths.get("temp/" + gene_type + "/lfcPair.txt");
        Calculator.createFileIfNotExists(path);

        try (BufferedWriter writer = Files.newBufferedWriter(path, StandardCharsets.UTF_8)) {
            StringBuilder stringBuilder = new StringBuilder();
            int pairCount = 0; // 计数器，用于控制写入频率

            // 使用普通循环创建基因对组合并写入文件
            for (int i = 0; i < uniqueGenesList.size(); i++) {
                for (int j = i; j < uniqueGenesList.size(); j++) {
                    // 生成基因对
                    stringBuilder.append(uniqueGenesList.get(i))
                            .append(":")
                            .append(uniqueGenesList.get(j))
                            .append(System.lineSeparator());
                    pairCount++;

                    // 每1000对写入一次
                    if (pairCount % 1000 == 0) {
                        writer.write(stringBuilder.toString());
                        stringBuilder.setLength(0); // 清空StringBuilder
                        writer.flush(); // 刷新缓冲区
                    }
                }
            }

            // 写入剩余的基因对
            if (stringBuilder.length() > 0) {
                writer.write(stringBuilder.toString());
                writer.flush(); // 刷新缓冲区
            }
        }
        return ecNumGeneHashMap;
    }
}
