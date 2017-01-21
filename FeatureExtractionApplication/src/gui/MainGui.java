package gui;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.EventQueue;
import java.awt.FlowLayout;
import java.awt.GridLayout;
import java.awt.Image;
import java.awt.TextArea;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;

import javax.imageio.ImageIO;
import javax.swing.DefaultComboBoxModel;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTabbedPane;
import javax.swing.ScrollPaneConstants;
import javax.swing.border.BevelBorder;

import org.opencv.features2d.KeyPoint;

import config.ImageData;
import imgproc.FeatureDescriptionImage;
import imgproc.FeatureExtractionImage2;
import imgproc.PanoramaImage;

public class MainGui {

	private JFrame frame;
	private JPanel panel_select_file;
	JLabel label_img1, label_img2, label_descriptorResult, label_panorama;

	/**
	 * Launch the application.
	 */
	public static void main(String[] args) {
		EventQueue.invokeLater(new Runnable() {
			public void run() {
				try {
					MainGui window = new MainGui();
					window.frame.setVisible(true);
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		});
	}

	/**
	 * Create the application.
	 */
	public MainGui() {
		initialize();
	}

	/**
	 * Initialize the contents of the frame.
	 */
	private void initialize() {
		frame = new JFrame();
		frame.setBounds(100, 100, 800, 550);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

		JTabbedPane tabbedPane = new JTabbedPane(JTabbedPane.TOP);
		frame.getContentPane().add(tabbedPane, BorderLayout.CENTER);

		JPanel panel_main = new JPanel();
		tabbedPane.addTab("Main", null, panel_main, null);
		panel_main.setLayout(new BorderLayout(0, 0));

		panel_select_file = new JPanel();
		panel_select_file.setSize(new Dimension(0, 100));
		FlowLayout fl_panel_select_file = (FlowLayout) panel_select_file.getLayout();
		panel_main.add(panel_select_file, BorderLayout.NORTH);

		JPanel panel_images = new JPanel();
		panel_main.add(panel_images, BorderLayout.CENTER);
		panel_images.setLayout(new GridLayout(0, 2, 0, 0));

		label_img1 = new JLabel("");
		label_img1.addMouseListener(new MouseAdapter() {
			@Override
			public void mouseClicked(MouseEvent e) {
				if (ImageData.feaExtImg1 != null && ImageData.feaExtImg2 != null) {
					StringBuilder sb = new StringBuilder();
					for (KeyPoint pt : ImageData.feaExtImg1.getKeyPoints().toArray()) {
						sb.append(pt + "\r\n");
					}
					showDialog(sb.toString());
				}
			}
		});
		label_img1.setBorder(new BevelBorder(BevelBorder.LOWERED, null, null, null, null));
		panel_images.add(label_img1);

		label_img2 = new JLabel("");
		label_img2.addMouseListener(new MouseAdapter() {
			@Override
			public void mouseClicked(MouseEvent e) {
				if (ImageData.feaExtImg1 != null && ImageData.feaExtImg2 != null) {
					StringBuilder sb = new StringBuilder();
					for (KeyPoint pt : ImageData.feaExtImg2.getKeyPoints().toArray()) {
						sb.append(pt + "\r\n");
					}
					showDialog(sb.toString());
				}
			}
		});
		label_img2.setBorder(new BevelBorder(BevelBorder.LOWERED, null, null, null, null));
		panel_images.add(label_img2);

		JPanel panel = new JPanel();
		panel_main.add(panel, BorderLayout.SOUTH);
		panel.setLayout(new FlowLayout(FlowLayout.CENTER, 5, 5));

		JButton button = new JButton("Original Image");
		button.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				showImage(ImageData.img1, label_img1);
				showImage(ImageData.img2, label_img2);
			}
		});
		panel.add(button);

		JButton button_1 = new JButton("1. Extract Feature");
		button_1.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				if (ImageData.img1 == null || ImageData.img2 == null) {
					JOptionPane.showMessageDialog(frame, "Please load image first");
					return;
				}

				// Feature Extraction on image
				ImageData.feaExtImg1 = new FeatureExtractionImage2(ImageData.img1.getAbsolutePath());
				ImageData.feaExtImg2 = new FeatureExtractionImage2(ImageData.img2.getAbsolutePath());

				// Show Image with KeyPoints
				showImage(ImageData.feaExtImg1.getBufferedImageKeypoints(), label_img1);
				showImage(ImageData.feaExtImg2.getBufferedImageKeypoints(), label_img2);
			}
		});
		panel.add(button_1);

		JButton btnDescribeFeature = new JButton("2. Describe Feature");
		btnDescribeFeature.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				if (ImageData.feaExtImg1 == null || ImageData.feaExtImg2 == null) {
					JOptionPane.showMessageDialog(frame, "Please extract feature first");
					return;
				}
				ImageData.feaDesImg = new FeatureDescriptionImage(ImageData.feaExtImg1, ImageData.feaExtImg2);
				showImage(ImageData.feaDesImg.getBufferedImgMatches(), label_descriptorResult);
				tabbedPane.setSelectedIndex(1);
			}
		});
		panel.add(btnDescribeFeature);

		JButton btnPanorama = new JButton("3. Panorama");
		btnPanorama.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				PanoramaImage img = new PanoramaImage(ImageData.feaExtImg1, ImageData.feaExtImg2, ImageData.feaDesImg);
				showImage(img.getBufferedPanoramaImage(), label_panorama);
				tabbedPane.setSelectedIndex(2);
			}
		});
		panel.add(btnPanorama);

		JPanel panel_descriptor = new JPanel();
		tabbedPane.addTab("Descriptor", null, panel_descriptor, null);
		panel_descriptor.setLayout(new BorderLayout(0, 0));

		label_descriptorResult = new JLabel("");
		panel_descriptor.add(label_descriptorResult);

		JPanel panel_panorama = new JPanel();
		tabbedPane.addTab("Panorama", null, panel_panorama, null);
		panel_panorama.setLayout(new BorderLayout(0, 0));

		label_panorama = new JLabel("");
		panel_panorama.add(label_panorama);

		JButton btnLoadImage_1 = new JButton("Load Image 1");
		btnLoadImage_1.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				ImageData.img1 = showFileSelector(e);
				showImage(ImageData.img1, label_img1);
			}
		});

		JComboBox comboBox = new JComboBox();
		comboBox.addItemListener(new ItemListener() {
			public void itemStateChanged(ItemEvent arg0) {

				if (!arg0.getItem().toString().contains("---")) {
					ImageData.img1 = new File(arg0.getItem().toString().split("/")[0]);
					ImageData.img2 = new File(arg0.getItem().toString().split("/")[1]);
					showImage(ImageData.img1, label_img1);
					showImage(ImageData.img2, label_img2);
				}
			}
		});
		comboBox.setModel(new DefaultComboBoxModel(new String[] { "---Please Select---", "po1.png/po2.png",
				"chu1.png/chu3.png", "till.png/box.png", "pan1.jpg/pan2.jpg", "OLAV.jpg/OLAV1.jpg", "p1.jpg/p2.jpg" }));
		comboBox.setToolTipText("");
		panel_select_file.add(comboBox);
		panel_select_file.add(btnLoadImage_1);

		JButton btnLoadImage_2 = new JButton("Load Image 2");
		btnLoadImage_2.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				ImageData.img2 = showFileSelector(e);
				showImage(ImageData.img2, label_img2);
			}
		});
		panel_select_file.add(btnLoadImage_2);
	}

	private File showFileSelector(ActionEvent e) {
		final JFileChooser fc = new JFileChooser();
		fc.setCurrentDirectory(new File(System.getProperty("user.dir")));
		int returnVal = fc.showOpenDialog((Component) e.getSource());
		if (returnVal == JFileChooser.APPROVE_OPTION) {
			return fc.getSelectedFile();
		}
		return null;
	}

	private void showImage(File imgFile, JLabel label) {

		if (imgFile == null)
			return;

		try {
			BufferedImage img = ImageIO.read(imgFile);
			showImage(img, label);

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private void showImage(BufferedImage img, JLabel label) {
		float labelRatio = (float) label.getWidth() / label.getHeight();
		float imgRatio = (float) img.getWidth() / img.getHeight();

		if (imgRatio > labelRatio)
			label.setIcon(new ImageIcon(
					img.getScaledInstance(label.getWidth(), (int) (label.getWidth() / imgRatio), Image.SCALE_SMOOTH)));
		else
			label.setIcon(new ImageIcon(img.getScaledInstance((int) (label.getHeight() * imgRatio), label.getHeight(),
					Image.SCALE_SMOOTH)));
	}

	private void showDialog(String message) {
		TextArea textArea = new TextArea();
		textArea.setText(message);

		JScrollPane scrollPane = new JScrollPane(textArea);
		scrollPane.setVerticalScrollBarPolicy(ScrollPaneConstants.VERTICAL_SCROLLBAR_ALWAYS);

		JOptionPane optionPane = new JOptionPane();
		optionPane.add(scrollPane);

		optionPane.showMessageDialog(null, textArea, "KeyPoints", JOptionPane.INFORMATION_MESSAGE);
	}

}
