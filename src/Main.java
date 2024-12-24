import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;


public class Main {
    public static void main(String[] args) throws Exception {
        // 0、初始化相关参数
        // String gene_type = "Human"; "Arabidopsis_thaliana"
        String gene_type = "Yeast";
        double alpha = 0.9D;

        // 1、随机游走的相关性得分矩阵
//        Calculator.getRWRMatrix(gene_type, alpha);
        System.out.println("矩阵已算完！");
        System.out.println("======================================");

//         2、计算基因相似度
//        LfcPair.getLfc(gene_type);
        String input = "temp/" + gene_type + "/lfcPair.txt";
        final String output = "result/" + gene_type + "/similarityResult" + alpha + ".txt";
        final Path outpath = Paths.get(output);
        Calculator.createFileIfNotExists(outpath);
        // 使用 TRUNCATE_EXISTING 确保覆盖原有内容
        Files.write(outpath, new byte[0], StandardOpenOption.TRUNCATE_EXISTING); // 先清空文件内容
        Calculator.data_initializer(gene_type, alpha);
        System.out.println("数据初始化完毕！");

        // 读取待计算的和已有结果的行
        String[] lines = Files.readAllLines(Paths.get(input)).toArray(new String[0]);
        int thCnt = Runtime.getRuntime().availableProcessors(); // 线程数
        ExecutorService executor = Executors.newFixedThreadPool(thCnt);
        CountDownLatch latch = new CountDownLatch(lines.length); // 计数器
        ConcurrentLinkedQueue<String> results = new ConcurrentLinkedQueue<>(); // 存储结果
        System.out.println("基因相似度计算中...");
        long start = System.nanoTime();
        for (String line : lines) {
            executor.submit(() -> {
                try {
                    String[] parts = line.split(":");
                    String g1 = parts[0];
                    String g2 = parts[1];
                    double sim = g1.equals(g2) ? 1.0 : Calculator.get_similarity_of_genes(g1, g2);
                    //若sim为NaN，则说明两基因没有共同的GO term，则直接输出0
                    if (Double.isNaN(sim)) {
                        sim = 0.0;
                        System.out.println(g1 + "和" + g2 + "没有共同的GO term，相似度为0.0");
                    }
                    // 收集结果
                    String result = g1 + "\t" + g2 + "\t" + sim;
                    results.add(result);
                    // 输出日志
                    if (results.size() % 10000 == 0) {
                        long end = System.nanoTime();
                        System.out.print(result + " 已处理:" + results.size() + "耗时:" + ((end - start) / 1000000) + "\n");
                    }
                } finally {
                    latch.countDown(); // 确保计数减少
                }
            });
        }
        // 等待所有任务完成
        latch.await();
        executor.shutdown(); // 先关闭 ExecutorService，再关闭线程池
        // 将结果写入文件
        try (BufferedWriter writer = Files.newBufferedWriter(outpath, StandardOpenOption.CREATE, StandardOpenOption.APPEND)) {
            while (!results.isEmpty()) {
                String result = results.poll(); // 从队列中取出结果
                if (result != null) {
                    writer.write(result);
                    writer.newLine();
                }
            }
        } catch (IOException e) {
            System.err.println("写入文件时发生异常: " + e.getMessage());
        }
        System.out.println("基因相似度计算完成！");

        // 3、根据ec分组计算各组的LFC得分
        Calculator.Evaluator.LFC(gene_type, alpha);
        System.out.println("LFC得分已算完！");
        System.out.println("Successfully!");
    }
}
