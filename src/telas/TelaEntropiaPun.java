/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package telas;

import java.awt.Desktop;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.imageio.ImageIO;
import javax.swing.ImageIcon;
import javax.swing.JFileChooser;
import javax.swing.JOptionPane;
import java.net.MalformedURLException;
import java.net.URL;
import javax.swing.filechooser.FileNameExtensionFilter;
import processamento.ProcessamentoImagem;


/**
 *
 * @author Pedro Henrique
 */
public class TelaEntropiaPun extends javax.swing.JFrame {
    private File pastaOrigem, pastaDestino;
    private boolean escolheuOrigem,escolheuDestino = false;
    private int numImgsProcessadas = 0;
    private int numImgsTotais;
    private boolean finalizado = false;
    /**
     * Creates new form TelaEntropiaPun
     */
    public TelaEntropiaPun(){
        initComponents();
//        ImageIcon icon = new ImageIcon("src/telas/TomatoIconBin_RedTomatoLaC.png"); 
//        icon.setImage(icon.getImage().getScaledInstance(tomatoIconBin.getWidth(), tomatoIconBin.getHeight(), 1));
//        tomatoIconBin.setIcon(icon);
//        tomatoIconBin2.setIcon(icon);
        setDefaultCloseOperation(javax.swing.WindowConstants.HIDE_ON_CLOSE);
        txtAvisos.setVisible(false);
    }

    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        btnSelPastaDestino = new javax.swing.JButton();
        jLabel4 = new javax.swing.JLabel();
        btnCancelar = new javax.swing.JButton();
        jLabel1 = new javax.swing.JLabel();
        txtAvisos = new javax.swing.JLabel();
        jLabel2 = new javax.swing.JLabel();
        txtPastaOrigem = new javax.swing.JTextField();
        btnSelPastaOrigem = new javax.swing.JButton();
        jLabel3 = new javax.swing.JLabel();
        txtPastaDestino = new javax.swing.JTextField();
        btnSelPastaDestino1 = new javax.swing.JButton();
        jLabel5 = new javax.swing.JLabel();
        btnPngOption = new javax.swing.JRadioButton();
        btnJpgOption = new javax.swing.JRadioButton();
        btnBmpOption = new javax.swing.JRadioButton();
        btnCancelar1 = new javax.swing.JButton();
        btnDetectar = new javax.swing.JButton();
        tomatoIconBin = new javax.swing.JLabel();
        tomatoIconBin2 = new javax.swing.JLabel();

        btnSelPastaDestino.setText("Selecionar...");
        btnSelPastaDestino.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnSelPastaDestinoActionPerformed(evt);
            }
        });

        jLabel4.setHorizontalAlignment(javax.swing.SwingConstants.CENTER);
        jLabel4.setText("Selecione as extensões das imagens que deseja realizar a deteccão:");

        btnCancelar.setText("Cancelar");
        btnCancelar.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnCancelarActionPerformed(evt);
            }
        });

        setDefaultCloseOperation(javax.swing.WindowConstants.EXIT_ON_CLOSE);

        jLabel1.setFont(new java.awt.Font("Ubuntu", 1, 20)); // NOI18N
        jLabel1.setHorizontalAlignment(javax.swing.SwingConstants.CENTER);
        jLabel1.setText("Entropia de Pun em Imagem de uma Pasta");

        txtAvisos.setForeground(new java.awt.Color(255, 0, 0));
        txtAvisos.setHorizontalAlignment(javax.swing.SwingConstants.CENTER);
        txtAvisos.setText("jLabel5");

        jLabel2.setHorizontalAlignment(javax.swing.SwingConstants.CENTER);
        jLabel2.setText("Selecione a pasta onde os arquivos de imagem estão contidos:");

        txtPastaOrigem.setText("Selecione uma pasta de origem...");

        btnSelPastaOrigem.setText("Selecionar...");
        btnSelPastaOrigem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnSelPastaOrigemActionPerformed(evt);
            }
        });

        jLabel3.setHorizontalAlignment(javax.swing.SwingConstants.CENTER);
        jLabel3.setText("Selecione a pasta onde deseja salvar as imagens processadas:");

        txtPastaDestino.setText("Selecione uma pasta de destino...");
        txtPastaDestino.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                txtPastaDestinoActionPerformed(evt);
            }
        });

        btnSelPastaDestino1.setText("Selecionar...");
        btnSelPastaDestino1.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnSelPastaDestino1ActionPerformed(evt);
            }
        });

        jLabel5.setHorizontalAlignment(javax.swing.SwingConstants.CENTER);
        jLabel5.setText("Selecione as extensões das imagens que deseja realizar a deteccão:");

        btnPngOption.setSelected(true);
        btnPngOption.setText("PNG");

        btnJpgOption.setSelected(true);
        btnJpgOption.setText("JPEG/JPG");

        btnBmpOption.setSelected(true);
        btnBmpOption.setText("BMP");

        btnCancelar1.setText("Cancelar");
        btnCancelar1.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnCancelar1ActionPerformed(evt);
            }
        });

        btnDetectar.setText("Detectar!");
        btnDetectar.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnDetectarActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(txtAvisos, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addContainerGap())
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                .addGap(0, 0, Short.MAX_VALUE)
                .addComponent(txtPastaDestino, javax.swing.GroupLayout.PREFERRED_SIZE, 468, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(btnSelPastaDestino1, javax.swing.GroupLayout.PREFERRED_SIZE, 128, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(45, 45, 45))
            .addGroup(layout.createSequentialGroup()
                .addGap(63, 63, 63)
                .addComponent(txtPastaOrigem, javax.swing.GroupLayout.PREFERRED_SIZE, 467, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addComponent(btnSelPastaOrigem, javax.swing.GroupLayout.PREFERRED_SIZE, 128, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(45, 45, 45))
            .addGroup(layout.createSequentialGroup()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createSequentialGroup()
                        .addGap(219, 219, 219)
                        .addComponent(btnPngOption)
                        .addGap(86, 86, 86)
                        .addComponent(btnJpgOption)
                        .addGap(88, 88, 88)
                        .addComponent(btnBmpOption))
                    .addGroup(layout.createSequentialGroup()
                        .addGap(56, 56, 56)
                        .addComponent(tomatoIconBin, javax.swing.GroupLayout.PREFERRED_SIZE, 54, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addGap(18, 18, 18)
                        .addComponent(jLabel1, javax.swing.GroupLayout.PREFERRED_SIZE, 458, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addGap(18, 18, 18)
                        .addComponent(tomatoIconBin2, javax.swing.GroupLayout.PREFERRED_SIZE, 54, javax.swing.GroupLayout.PREFERRED_SIZE)))
                .addContainerGap(56, Short.MAX_VALUE))
            .addGroup(layout.createSequentialGroup()
                .addGap(110, 110, 110)
                .addComponent(btnCancelar1, javax.swing.GroupLayout.PREFERRED_SIZE, 160, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addComponent(btnDetectar, javax.swing.GroupLayout.PREFERRED_SIZE, 160, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(94, 94, 94))
            .addGroup(layout.createSequentialGroup()
                .addGap(20, 20, 20)
                .addComponent(jLabel3, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addContainerGap())
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(jLabel2, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addContainerGap())
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(jLabel5, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addContainerGap())
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addGap(18, 18, 18)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(tomatoIconBin2, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.PREFERRED_SIZE, 52, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addGroup(layout.createSequentialGroup()
                        .addGap(9, 9, 9)
                        .addComponent(jLabel1))
                    .addComponent(tomatoIconBin, javax.swing.GroupLayout.PREFERRED_SIZE, 52, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(jLabel2)
                .addGap(18, 18, 18)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(txtPastaOrigem, javax.swing.GroupLayout.PREFERRED_SIZE, 30, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(btnSelPastaOrigem, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                .addGap(18, 18, 18)
                .addComponent(jLabel3)
                .addGap(18, 18, 18)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(txtPastaDestino, javax.swing.GroupLayout.PREFERRED_SIZE, 30, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(btnSelPastaDestino1, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                .addGap(18, 18, 18)
                .addComponent(jLabel5)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(btnPngOption)
                    .addComponent(btnJpgOption)
                    .addComponent(btnBmpOption))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(txtAvisos)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(btnCancelar1)
                    .addComponent(btnDetectar))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        pack();
    }// </editor-fold>//GEN-END:initComponents

    private void btnSelPastaOrigemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnSelPastaOrigemActionPerformed
        // TODO add your handling code here:
        try{
            JFileChooser fs = new JFileChooser(new File("/home/"));
            fs.setDialogTitle("Selecione o Arquivo de Imagem");
            fs.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
            fs.setAcceptAllFileFilterUsed(false);
            int retorno = fs.showOpenDialog(null);
            if (retorno == JFileChooser.APPROVE_OPTION){
                txtPastaOrigem.setText(fs.getSelectedFile().getPath());
                pastaOrigem = fs.getSelectedFile();
                escolheuOrigem = true;
            }
        } catch(NullPointerException ex){
            Logger.getLogger(TelaDetectaPelePasta.class.getName()).log(Level.SEVERE, null, ex);
        }
    }//GEN-LAST:event_btnSelPastaOrigemActionPerformed

    private void txtPastaDestinoActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_txtPastaDestinoActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_txtPastaDestinoActionPerformed

    private void btnSelPastaDestinoActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnSelPastaDestinoActionPerformed
        // TODO add your handling code here:
        try{
            JFileChooser fs = new JFileChooser(new File("/home/"));
            fs.setDialogTitle("Selecione o Arquivo de Imagem");
            fs.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
            fs.setAcceptAllFileFilterUsed(false);
            int retorno = fs.showOpenDialog(null);
            if (retorno == JFileChooser.APPROVE_OPTION){
                txtPastaDestino.setText(fs.getSelectedFile().getPath());
                pastaDestino = fs.getSelectedFile();
                escolheuDestino = true;
            }
        } catch(NullPointerException ex){
            Logger.getLogger(TelaDetectaPelePasta.class.getName()).log(Level.SEVERE, null, ex);
        }
    }//GEN-LAST:event_btnSelPastaDestinoActionPerformed

    private void btnSelPastaDestino1ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnSelPastaDestino1ActionPerformed
        // TODO add your handling code here:
        try{
            JFileChooser fs = new JFileChooser(new File("/home/"));
            fs.setDialogTitle("Selecione o Arquivo de Imagem");
            fs.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
            fs.setAcceptAllFileFilterUsed(false);
            int retorno = fs.showOpenDialog(null);
            if (retorno == JFileChooser.APPROVE_OPTION){
                txtPastaDestino.setText(fs.getSelectedFile().getPath());
                pastaDestino = fs.getSelectedFile();
                escolheuDestino = true;
            }
        } catch(NullPointerException ex){
            Logger.getLogger(TelaDetectaPelePasta.class.getName()).log(Level.SEVERE, null, ex);
        }
    }//GEN-LAST:event_btnSelPastaDestino1ActionPerformed

    private void btnCancelarActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnCancelarActionPerformed
        // TODO add your handling code here:
        this.dispose();
    }//GEN-LAST:event_btnCancelarActionPerformed

    private void btnCancelar1ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnCancelar1ActionPerformed
        // TODO add your handling code here:
        this.dispose();
    }//GEN-LAST:event_btnCancelar1ActionPerformed

    private void btnDetectarActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnDetectarActionPerformed
        // TODO add your handling code here:
        if(escolheuOrigem && escolheuDestino){
            finalizado = false;
            numImgsProcessadas = 0;
            numImgsTotais = 0;
            txtAvisos.setVisible(false);
            for(File arqPastaOrigem: pastaOrigem.listFiles()){
                if((btnPngOption.isSelected() && arqPastaOrigem.getName().toLowerCase().endsWith(".png"))||
                    (btnJpgOption.isSelected() && arqPastaOrigem.getName().toLowerCase().endsWith(".jpg")) ||
                    (btnJpgOption.isSelected() && arqPastaOrigem.getName().toLowerCase().endsWith(".jpeg")) ||
                    (btnBmpOption.isSelected() && arqPastaOrigem.getName().toLowerCase().endsWith(".bmp")))
                {
                    numImgsTotais++;
                }
            }
            Thread t1 = new Thread(){
                @Override
                public void run()
                {
                    for(File arqPastaOrigem: pastaOrigem.listFiles()){
                        if((btnPngOption.isSelected() && arqPastaOrigem.getName().toLowerCase().endsWith(".png"))||
                            (btnJpgOption.isSelected() && arqPastaOrigem.getName().toLowerCase().endsWith(".jpg")) ||
                            (btnJpgOption.isSelected() && arqPastaOrigem.getName().toLowerCase().endsWith(".jpeg")) ||
                            (btnBmpOption.isSelected() && arqPastaOrigem.getName().toLowerCase().endsWith(".bmp")))
                        {
                            try {
                                BufferedImage imgOriginal = ImageIO.read(arqPastaOrigem);
                                BufferedImage imagemProcessada = ProcessamentoImagem.EntropiaPun(imgOriginal);
                                
                                File fo;
                                if(arqPastaOrigem.getName().length() >= 16){
                                    String finalRedness = "_RedCIRG2.png";
                                    if(finalRedness.equals(arqPastaOrigem.getName().substring(arqPastaOrigem.getName().length()-16)))
                                        fo = new File(pastaDestino.getPath()+"/"+arqPastaOrigem.getName().substring(0, arqPastaOrigem.getName().length() - 17)+"_RedCIRG2_BinPun"+".png");
                                    else
                                        fo = new File(pastaDestino.getPath()+"/"+arqPastaOrigem.getName().substring(0, arqPastaOrigem.getName().length() - 4)+"_RedCIRG2_BinPun"+".png");
                                }
                                else{
                                    fo = new File(pastaDestino.getPath()+"/"+arqPastaOrigem.getName().substring(0, arqPastaOrigem.getName().length() - 4)+"_RedCIRG2_BinPun"+".png");
                                }
                                ImageIO.write(imagemProcessada, "png", fo);
                                numImgsProcessadas++;
                            } catch (IOException ex) {
                                Logger.getLogger(TelaDetectaPelePasta.class.getName()).log(Level.SEVERE, null, ex);
                            }
                        }
                    }
                    finalizado = true;
                }
            };
            t1.start();

            Thread t2 = new Thread(){
                @Override
                public void run()
                {
                    TelaCarregamento tc = new TelaCarregamento();
                    tc.setVisible(true);
                    while(!finalizado){
                        tc.AtualizaCarregamento(numImgsProcessadas, numImgsTotais);
                    }
                    tc.dispose();
                    if(numImgsProcessadas > 0 ){
                        JOptionPane.showMessageDialog(rootPane, numImgsProcessadas+" imagens foram processadas com sucesso!", "Detecção Realizada!", JOptionPane.INFORMATION_MESSAGE);
                        Desktop desktop = Desktop.getDesktop();
                        try {
                            desktop.open(pastaDestino);
                        }catch (IOException ex) {
                            Logger.getLogger(TelaDetectaPelePasta.class.getName()).log(Level.SEVERE, null, ex);
                        }
                    }else{
                        txtAvisos.setText("Não foram encontradas imagens, com as extensões especificadas, na pasta de origem!");
                        txtAvisos.setVisible(true);
                    }
                }
            };
            t2.start();

        }else if(!escolheuOrigem){
            txtAvisos.setText("A pasta de origem não foi especificada!");
            txtAvisos.setVisible(true);
        }else if(!escolheuDestino){
            txtAvisos.setText("A pasta de destino não foi especificada!");
            txtAvisos.setVisible(true);
        }
    }//GEN-LAST:event_btnDetectarActionPerformed

    /**
     * @param args the command line arguments
     */
    public static void main(String args[]) throws MalformedURLException {
        /* Set the Nimbus look and feel */
        //<editor-fold defaultstate="collapsed" desc=" Look and feel setting code (optional) ">
        /* If Nimbus (introduced in Java SE 6) is not available, stay with the default look and feel.
         * For details see http://download.oracle.com/javase/tutorial/uiswing/lookandfeel/plaf.html 
         */
        try {
            for (javax.swing.UIManager.LookAndFeelInfo info : javax.swing.UIManager.getInstalledLookAndFeels()) {
                if ("Nimbus".equals(info.getName())) {
                    javax.swing.UIManager.setLookAndFeel(info.getClassName());
                    break;
                }
            }
        } catch (ClassNotFoundException ex) {
            java.util.logging.Logger.getLogger(TelaEntropiaPun.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (InstantiationException ex) {
            java.util.logging.Logger.getLogger(TelaEntropiaPun.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (IllegalAccessException ex) {
            java.util.logging.Logger.getLogger(TelaEntropiaPun.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (javax.swing.UnsupportedLookAndFeelException ex) {
            java.util.logging.Logger.getLogger(TelaEntropiaPun.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        }
        //</editor-fold>
        //</editor-fold>

        /* Create and display the form */
        java.awt.EventQueue.invokeLater(new Runnable() {
            public void run() {
                new TelaEntropiaPun().setVisible(true);
            }
        });
    }

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JRadioButton btnBmpOption;
    private javax.swing.JButton btnCancelar;
    private javax.swing.JButton btnCancelar1;
    private javax.swing.JButton btnDetectar;
    private javax.swing.JRadioButton btnJpgOption;
    private javax.swing.JRadioButton btnPngOption;
    private javax.swing.JButton btnSelPastaDestino;
    private javax.swing.JButton btnSelPastaDestino1;
    private javax.swing.JButton btnSelPastaOrigem;
    private javax.swing.JLabel jLabel1;
    private javax.swing.JLabel jLabel2;
    private javax.swing.JLabel jLabel3;
    private javax.swing.JLabel jLabel4;
    private javax.swing.JLabel jLabel5;
    private javax.swing.JLabel tomatoIconBin;
    private javax.swing.JLabel tomatoIconBin2;
    private javax.swing.JLabel txtAvisos;
    private javax.swing.JTextField txtPastaDestino;
    private javax.swing.JTextField txtPastaOrigem;
    // End of variables declaration//GEN-END:variables
}
