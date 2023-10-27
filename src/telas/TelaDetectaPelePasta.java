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
import javax.swing.JFileChooser;
import javax.swing.JOptionPane;
import java.net.MalformedURLException;
import java.net.URL;
import javax.swing.ImageIcon;
import javax.swing.filechooser.FileNameExtensionFilter;
import processamento.ProcessamentoImagem;

/**
 *
 * @author Álef e Ketson
 */
public class TelaDetectaPelePasta extends javax.swing.JFrame {
    private File pastaOrigem, pastaDestino;
    private boolean escolheuOrigem,escolheuDestino = false;
    private int numImgsProcessadas = 0;
    private int numImgsTotais;
    private boolean finalizado = false;

    /**
     * Creates new form TelaPostPasta
     */
    public TelaDetectaPelePasta() {        
        initComponents();
//        ImageIcon icon = new ImageIcon("src/telas/TomatoIcon.png"); ;
//        icon.setImage(icon.getImage().getScaledInstance(tomatoIcon.getWidth(), tomatoIcon.getHeight(), 1));
//        tomatoIcon.setIcon(icon);
//        tomatoIcon2.setIcon(icon);
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

        jPanel1 = new javax.swing.JPanel();
        jLabel1 = new javax.swing.JLabel();
        jLabel2 = new javax.swing.JLabel();
        btnSelPastaOrigem = new javax.swing.JButton();
        jLabel3 = new javax.swing.JLabel();
        btnSelPastaDestino = new javax.swing.JButton();
        jLabel4 = new javax.swing.JLabel();
        btnPngOption = new javax.swing.JRadioButton();
        btnJpgOption = new javax.swing.JRadioButton();
        btnBmpOption = new javax.swing.JRadioButton();
        btnDetectar = new javax.swing.JButton();
        txtAvisos = new javax.swing.JLabel();
        btnCancelar = new javax.swing.JButton();
        txtPastaOrigem = new javax.swing.JTextField();
        txtPastaDestino = new javax.swing.JTextField();

        javax.swing.GroupLayout jPanel1Layout = new javax.swing.GroupLayout(jPanel1);
        jPanel1.setLayout(jPanel1Layout);
        jPanel1Layout.setHorizontalGroup(
            jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 100, Short.MAX_VALUE)
        );
        jPanel1Layout.setVerticalGroup(
            jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 100, Short.MAX_VALUE)
        );

        setDefaultCloseOperation(javax.swing.WindowConstants.DISPOSE_ON_CLOSE);
        setFocusCycleRoot(false);
        setResizable(false);

        jLabel1.setFont(new java.awt.Font("Ubuntu", 1, 20)); // NOI18N
        jLabel1.setForeground(new java.awt.Color(255, 44, 6));
        jLabel1.setHorizontalAlignment(javax.swing.SwingConstants.CENTER);
        jLabel1.setText("Detecção de Redness em Imagens de uma Pasta");

        jLabel2.setHorizontalAlignment(javax.swing.SwingConstants.CENTER);
        jLabel2.setText("Selecione a pasta onde os arquivos de imagem estão contidos:");

        btnSelPastaOrigem.setText("Selecionar...");
        btnSelPastaOrigem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnSelPastaOrigemActionPerformed(evt);
            }
        });

        jLabel3.setHorizontalAlignment(javax.swing.SwingConstants.CENTER);
        jLabel3.setText("Selecione a pasta onde deseja salvar as imagens processadas:");

        btnSelPastaDestino.setText("Selecionar...");
        btnSelPastaDestino.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnSelPastaDestinoActionPerformed(evt);
            }
        });

        jLabel4.setHorizontalAlignment(javax.swing.SwingConstants.CENTER);
        jLabel4.setText("Selecione as extensões das imagens que deseja realizar a deteccão:");

        btnPngOption.setSelected(true);
        btnPngOption.setText("PNG");

        btnJpgOption.setSelected(true);
        btnJpgOption.setText("JPEG/JPG");

        btnBmpOption.setSelected(true);
        btnBmpOption.setText("BMP");

        btnDetectar.setText("Detectar!");
        btnDetectar.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnDetectarActionPerformed(evt);
            }
        });

        txtAvisos.setForeground(new java.awt.Color(255, 0, 0));
        txtAvisos.setHorizontalAlignment(javax.swing.SwingConstants.CENTER);
        txtAvisos.setText("jLabel5");

        btnCancelar.setText("Cancelar");
        btnCancelar.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnCancelarActionPerformed(evt);
            }
        });

        txtPastaOrigem.setText("Selecione uma pasta de origem...");

        txtPastaDestino.setText("Selecione uma pasta de origem...");
        txtPastaDestino.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                txtPastaDestinoActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addGap(55, 55, 55)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                        .addComponent(txtPastaOrigem)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(btnSelPastaOrigem, javax.swing.GroupLayout.PREFERRED_SIZE, 128, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addGap(54, 54, 54))
                    .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                        .addComponent(txtPastaDestino)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(btnSelPastaDestino, javax.swing.GroupLayout.PREFERRED_SIZE, 128, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addGap(53, 53, 53))))
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(jLabel3, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(jLabel4, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addGroup(layout.createSequentialGroup()
                        .addGap(0, 0, Short.MAX_VALUE)
                        .addComponent(jLabel1, javax.swing.GroupLayout.PREFERRED_SIZE, 474, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addGap(72, 72, 72)))
                .addContainerGap())
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                .addGap(0, 35, Short.MAX_VALUE)
                .addComponent(jLabel2, javax.swing.GroupLayout.PREFERRED_SIZE, 549, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(58, 58, 58))
            .addGroup(layout.createSequentialGroup()
                .addGap(155, 155, 155)
                .addComponent(btnPngOption)
                .addGap(84, 84, 84)
                .addComponent(btnJpgOption)
                .addGap(82, 82, 82)
                .addComponent(btnBmpOption)
                .addGap(0, 0, Short.MAX_VALUE))
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(txtAvisos, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addContainerGap())
            .addGroup(layout.createSequentialGroup()
                .addGap(84, 84, 84)
                .addComponent(btnCancelar, javax.swing.GroupLayout.PREFERRED_SIZE, 160, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addComponent(btnDetectar, javax.swing.GroupLayout.PREFERRED_SIZE, 160, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(80, 80, 80))
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addGap(30, 30, 30)
                .addComponent(jLabel1)
                .addGap(16, 16, 16)
                .addComponent(jLabel2)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(txtPastaOrigem, javax.swing.GroupLayout.PREFERRED_SIZE, 30, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(btnSelPastaOrigem, javax.swing.GroupLayout.PREFERRED_SIZE, 30, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(18, 18, 18)
                .addComponent(jLabel3)
                .addGap(18, 18, 18)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(txtPastaDestino, javax.swing.GroupLayout.PREFERRED_SIZE, 30, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(btnSelPastaDestino, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(jLabel4)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(btnPngOption)
                    .addComponent(btnJpgOption)
                    .addComponent(btnBmpOption))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addComponent(txtAvisos)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(btnCancelar)
                    .addComponent(btnDetectar))
                .addGap(15, 15, 15))
        );

        pack();
        setLocationRelativeTo(null);
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
                                BufferedImage imagemProcessada = ProcessamentoImagem.detectarPeleHumana(imgOriginal);
                                File fo = new File(pastaDestino.getPath()+"/"+arqPastaOrigem.getName().substring(0, arqPastaOrigem.getName().length() - 4)+"_RedCIRG2"+".png");
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

    private void btnCancelarActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnCancelarActionPerformed
        // TODO add your handling code here:
        this.dispose();
    }//GEN-LAST:event_btnCancelarActionPerformed

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

    private void txtPastaDestinoActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_txtPastaDestinoActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_txtPastaDestinoActionPerformed

    /**
     * @param args the command line arguments
     */
    public static void main(String args[]) {
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
            java.util.logging.Logger.getLogger(TelaDetectaPelePasta.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (InstantiationException ex) {
            java.util.logging.Logger.getLogger(TelaDetectaPelePasta.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (IllegalAccessException ex) {
            java.util.logging.Logger.getLogger(TelaDetectaPelePasta.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (javax.swing.UnsupportedLookAndFeelException ex) {
            java.util.logging.Logger.getLogger(TelaDetectaPelePasta.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        }
        //</editor-fold>
        //</editor-fold>
        //</editor-fold>
        //</editor-fold>
        //</editor-fold>
        //</editor-fold>
        //</editor-fold>
        //</editor-fold>

        /* Create and display the form */
        java.awt.EventQueue.invokeLater(new Runnable() {
            public void run() {
                new TelaDetectaPelePasta().setVisible(true);
            }
        });
    }

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JRadioButton btnBmpOption;
    private javax.swing.JButton btnCancelar;
    private javax.swing.JButton btnDetectar;
    private javax.swing.JRadioButton btnJpgOption;
    private javax.swing.JRadioButton btnPngOption;
    private javax.swing.JButton btnSelPastaDestino;
    private javax.swing.JButton btnSelPastaOrigem;
    private javax.swing.JLabel jLabel1;
    private javax.swing.JLabel jLabel2;
    private javax.swing.JLabel jLabel3;
    private javax.swing.JLabel jLabel4;
    private javax.swing.JPanel jPanel1;
    private javax.swing.JLabel txtAvisos;
    private javax.swing.JTextField txtPastaDestino;
    private javax.swing.JTextField txtPastaOrigem;
    // End of variables declaration//GEN-END:variables
}
