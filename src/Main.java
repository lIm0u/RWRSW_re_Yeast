import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;


public class Main {
    public static void main(String[] args) throws Exception {
        // 0、初始化相关参数
        // String gene_type = "Human"; "Arabidopsis_thaliana"
        String gene_type = "Yeast";
//        String gene_type = args[0];
        double alpha = 0.9D;

        // 1、随机游走的相关性得分矩阵
//        Calculator.getRWRMatrix(gene_type, alpha);
        System.out.println("矩阵已算完！");
        System.out.println("======================================");

        // 2、计算基因相似度
        String input = "buf/" + gene_type + "/lfc.txt";
        final String output = "result/" + gene_type + "/similarityResult" + alpha + ".txt";
        final Path outpath = Paths.get(output);
        Calculator.createFileIfNotExists(outpath);
        Calculator.data_initializer(gene_type, alpha);
        System.out.println("数据初始化完毕！");

        // 读取待计算的和已有结果的行
        String[] lines = Files.readAllLines(Paths.get(input)).toArray(new String[0]);
        final String[] lines2 = Files.readAllLines(outpath).toArray(new String[0]);
        // 若结果文件已有部分结果，则跳过已有结果
        if (lines2.length > 0) {
            // 使用 Set 存储已有结果以提高查找速度
            Set<String> existingResults = new HashSet<>();
            for (String line : lines2) {
                String[] parts = line.split("\t");
                if (parts.length >= 2) {
                    // 仅使用 g1 和 g2 生成组合键
                    String key = parts[0].trim() + ":" + parts[1].trim(); // 组合键
                    existingResults.add(key);
                }
            }
            // 过滤 lines，只保留那些不在 existingResults 中的行
            lines = Arrays.stream(lines)
                    .filter(line -> !existingResults.contains(line.trim()))
                    .toArray(String[]::new);
        }
        int thCnt = Runtime.getRuntime().availableProcessors(); // 线程数
        ExecutorService executor = Executors.newFixedThreadPool(thCnt);
        CountDownLatch latch = new CountDownLatch(lines.length); // 计数器，等待所有线程完成
        ConcurrentLinkedQueue<String> results = new ConcurrentLinkedQueue<>(); // 存储结果
        System.out.println("基因相似度计算中...");
        for (int i = 0; i < thCnt; i++) {
            final int startIndex = i;
            String[] finalLines = lines;
            executor.submit(() -> {
                try {
                    for (int j = startIndex; j < finalLines.length; j += thCnt) {
                        long taskStartTime = System.currentTimeMillis(); // 记录每个任务开始时间
                        String[] parts = finalLines[j].split(":");
                        if (parts.length < 2) {
                            continue; // 跳过当前无效行
                        }
                        String g1 = parts[0];
                        String g2 = parts[1];
                        double sim;
                        if (g1.equals(g2)) {
                            sim = 1.0;
                        } else {
                            sim = Calculator.get_similarity_of_genes(g1, g2);
                        }
                        // 收集结果
                        String result = g1 + "\t" + g2 + "\t" + sim;
                        long taskEndTime = System.currentTimeMillis(); // 记录每个任务结束时间
                        results.add(result);
                        System.out.print("Thread:" + (startIndex + 1) + "\t" + result + " Task:" + (j + 1) + "\tTime Use:" + (taskEndTime - taskStartTime) + "\n");
                    }
                } finally {
                    latch.countDown(); // 完成后减少计数
                }
            });
        }
        // 等待所有任务完成
        latch.await();
        executor.shutdown(); // 先关闭 ExecutorService，再关闭线程池
        // 将结果写入文件
        synchronized (outpath) {
            try (BufferedWriter writer = Files.newBufferedWriter(outpath, StandardOpenOption.CREATE, StandardOpenOption.APPEND)) {
                while (!results.isEmpty()) {
                    String result = results.poll(); // 从队列中取出结果
                    if (result != null) {
                        writer.write(result);
                        writer.newLine();
                    }
                }
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
        }
        System.out.println("基因相似度计算完成！");

//        // 3、根据ec分组计算各组的LFC得分
//        Calculator.Evaluator.LFC(gene_type, alpha);
//        System.out.println("LFC得分已算完！");
//        System.out.println("Successfully!");
    }
}
