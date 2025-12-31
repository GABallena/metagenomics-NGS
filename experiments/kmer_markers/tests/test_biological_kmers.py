import unittest
import numpy as np
import os
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import SGDClassifier
from sklearn.metrics import classification_report
from biological_kmers import (
    extract_features,
    create_feature_vector,
    generate_random_sequences,
    get_entropy,
    get_codon_bias,
    get_nucleotide_freq,
    get_repeat_content,
    get_kmer_entropy,
    get_palindrome_content,
    get_frame_bias,
    train_kmer_classifier,
    predict_kmer,
    parse_repeatmasker_output,
    parse_augustus_output,
    parse_trnascan_output,
    analyze_phylogenetic_signal,
    ModelPersistence,
    feature_generator
)
import itertools

# UniProt SwissProt database file
UNIPROT_SPROT_PATH = os.path.join(os.path.dirname(__file__), "data", "uniprot_sprot.fasta")

class TestBiologicalKmers(unittest.TestCase):
    # Standard kmer sizes for biological sequence analysis
    KMER_21 = 21
    KMER_31 = 31
    KMER_41 = 41
    KMER_51 = 51
    VALID_KMER_SIZES = [KMER_21, KMER_31, KMER_41, KMER_51]

    def setUp(self):
        # Ensure uniprot directory exists
        os.makedirs(os.path.dirname(UNIPROT_SPROT_PATH), exist_ok=True)
        # Test sequences
        self.test_seq = "ATGCATGCATGC"
        self.simple_seq = "AAAA"
        self.complex_seq = "ATGCTAGCTAGCTAG"
        self.palindrome_seq = "ATGCGCAT"

    def test_extract_features(self):
        features = extract_features(self.test_seq)
        self.assertEqual(len(features), 2)  # GC content and complexity
        self.assertAlmostEqual(features[0], 0.5)  # GC content
        self.assertAlmostEqual(features[1], 0.333, places=3)  # Complexity

    def test_create_feature_vector(self):
        # Test standard case
        vector = create_feature_vector(self.test_seq)
        self.assertIsInstance(vector, list)
        self.assertTrue(all(isinstance(x, (int, float)) for x in vector))
        
        # Test with known sequence
        known_seq = "ATGCATGC"
        vector = create_feature_vector(known_seq)
        self.assertEqual(vector[0], 0.5)  # GC content should be 0.5
        
        # Test with empty sequence
        with self.assertRaises(ValueError):
            create_feature_vector("")

    def test_generate_random_sequences(self):
        # Test with smallest standard kmer size
        seqs = generate_random_sequences(kmer_size=self.KMER_21)
        self.assertTrue(all(len(s) == self.KMER_21 for s in seqs))

        # Test with standard kmer sizes
        for k in self.VALID_KMER_SIZES[:2]:  # Only test first two to avoid memory issues
            seqs = generate_random_sequences(kmer_size=k)
            self.assertTrue(all(len(s) == k for s in seqs))

        # Test invalid kmer sizes
        with self.assertRaises(ValueError):
            generate_random_sequences(kmer_size=20)  # Not in valid sizes
        
        # Test warning for large kmer sizes
        with self.assertWarns(Warning):
            generate_random_sequences(kmer_size=self.KMER_41)

    def test_get_entropy(self):
        # Test uniform sequence
        self.assertEqual(get_entropy(self.simple_seq), 0.0)
        # Test diverse sequence
        self.assertGreater(get_entropy(self.complex_seq), 0.0)

    def test_get_codon_bias(self):
        # Test sequence shorter than 3
        self.assertEqual(get_codon_bias("AT"), 0)
        # Test normal sequence
        self.assertLessEqual(get_codon_bias(self.test_seq), 1.0)
        self.assertGreaterEqual(get_codon_bias(self.test_seq), 0.0)

    def test_get_nucleotide_freq(self):
        # Test dinucleotide frequencies
        freq = get_nucleotide_freq(self.test_seq, 2)
        self.assertIsInstance(freq, dict)
        self.assertEqual(sum(freq.values()), 1.0)

    def test_get_repeat_content(self):
        # Test sequence with repeats
        repeat_seq = "ATGATGATG"
        self.assertGreater(get_repeat_content(repeat_seq), 0.0)
        # Test sequence without repeats
        self.assertEqual(get_repeat_content(self.simple_seq), 1.0)

    def test_get_kmer_entropy(self):
        # Test with k=3
        entropy = get_kmer_entropy(self.test_seq, k=3)
        self.assertGreaterEqual(entropy, 0.0)
        # Test sequence shorter than k
        self.assertEqual(get_kmer_entropy("AT", k=3), 0.0)

    def test_get_palindrome_content(self):
        # Test palindrome sequence
        self.assertGreater(get_palindrome_content(self.palindrome_seq), 0.0)
        # Test non-palindrome sequence
        self.assertEqual(get_palindrome_content(self.simple_seq), 0.0)

    def test_get_frame_bias(self):
        # Test sequence with frame bias
        self.assertLessEqual(get_frame_bias(self.test_seq), 1.0)
        self.assertGreaterEqual(get_frame_bias(self.test_seq), 0.0)
        # Test short sequence
        self.assertEqual(get_frame_bias("AT"), 0)

    def test_edge_cases(self):
        # Test empty sequence
        self.assertEqual(get_entropy(""), 0.0)
        self.assertEqual(get_codon_bias(""), 0.0)
        self.assertEqual(len(generate_random_sequences(0)), 0)

    def test_input_validation(self):
        # Test invalid input
        with self.assertRaises(ValueError):
            extract_features("")
        with self.assertRaises(ValueError):
            get_nucleotide_freq(self.test_seq, 0)

    def test_train_kmer_classifier(self):
        try:
            # Train classifier
            clf, model_path = train_kmer_classifier()
            
            # Validate model structure and metrics
            self.assertIsInstance(clf, SGDClassifier)
            self.assertTrue(hasattr(clf, 'coef_'))
            self.assertTrue(isinstance(model_path, str))
            
            # Test incremental learning
            test_features = create_feature_vector(self.test_seq)
            prediction = clf.predict([test_features])
            self.assertIn(prediction[0], [0, 1])
            
        except ImportError as e:
            self.skipTest(f"Skipping due to missing dependencies: {str(e)}")
        except Exception as e:
            self.fail(f"Test failed with error: {str(e)}")

    def test_predict_kmer(self):
        try:
            # Test valid k-mer using standard size
            valid_kmer = "A" * self.KMER_31
            result = predict_kmer(valid_kmer)
            self.assertIn(result, ["Biological", "Artifact"])
            
            # Test with different valid kmer sizes
            for size in self.VALID_KMER_SIZES[:2]:  # Test first two sizes
                test_kmer = "A" * size
                result = predict_kmer(test_kmer)
                self.assertIn(result, ["Biological", "Artifact"])
                
            # Test invalid k-mer length
            with self.assertRaises(ValueError):
                predict_kmer("AT")  # Too short
                
            # Test invalid characters
            with self.assertRaises(ValueError):
                predict_kmer("N" * self.KMER_31)
                
        except Exception as e:
            self.fail(f"Test failed with error: {str(e)}")

    def test_external_tool_validation(self):
        # Test RepeatMasker output parsing
        repeat_output = """
        SW    perc perc perc  query     position in query    matching  repeat      position in repeat
        score  div. del. ins.  sequence  begin  end          repeat    class/family  begin  end
        
        1234   11.5  6.2  0.0  Seq1     100    200         AluY      SINE/Alu      1    100
        """
        features = parse_repeatmasker_output(repeat_output)
        self.assertEqual(features['repeat_count'], 1)
        self.assertEqual(len(features['repeat_types']), 1)
        self.assertEqual(features['repeat_lengths'][0], 100)
        
        # Test AUGUSTUS output parsing
        augustus_output = """
        # start gene g1
        Seq1    AUGUSTUS    gene    1    1000    1    +    .    g1
        Seq1    AUGUSTUS    CDS     1    500     1    +    0    g1
        """
        features = parse_augustus_output(augustus_output)
        self.assertEqual(features['gene_count'], 1)
        self.assertEqual(len(features['cds_lengths']), 1)
        self.assertEqual(features['cds_lengths'][0], 499)
        
        # Test tRNAscan-SE output parsing
        trna_output = """
        Sequence    tRNA    Bounds  tRNA    Anti    Intron Bounds   Inf
        Name        #      Begin    End     Type    Codon   Begin    End     Score
        --------    ----  ----    ----    ----    -----   ----    ----    -----
        Seq1        1     1000    1072    Leu     CAA     0        0       87.6
        """
        features = parse_trnascan_output(trna_output)
        self.assertEqual(features['trna_count'], 1)
        self.assertEqual(len(features['anticodon_types']), 1)
        self.assertEqual(features['scores'][0], 87.6)

    def test_phylogenetic_analysis(self):
        test_seq = "ATGCATGCATGC"
        signal = analyze_phylogenetic_signal(test_seq)
        self.assertGreaterEqual(signal, 0.0)
        self.assertLessEqual(signal, 1.0)
        
        # Test with known conserved sequence
        conserved_seq = "ATGGCCAAGTAA"  # Common start-stop pattern
        signal = analyze_phylogenetic_signal(conserved_seq)
        self.assertGreater(signal, 0.5)

    def test_model_persistence(self):
        try:
            # Create and fit a test model using SGD instead of RandomForest
            clf = SGDClassifier(random_state=42, loss='log_loss')  # Changed to SGD
            X = np.array([[0.5, 0.5], [0.3, 0.7], [0.8, 0.2]])
            y = np.array([1, 1, 0])
            clf.fit(X, y)
            
            # Test model saving and loading
            model_handler = ModelPersistence(model_dir="test_models")
            metadata = {'test': True}
            model_path = model_handler.save_model(clf, metadata)
            
            # Load and verify model
            loaded_clf = model_handler.load_model(model_path)
            np.testing.assert_array_equal(
                loaded_clf.predict([[0.5, 0.5]]),
                clf.predict([[0.5, 0.5]])
            )
            
            # Test loading latest model
            latest_clf = model_handler.load_latest_model()
            np.testing.assert_array_equal(
                latest_clf.predict([[0.5, 0.5]]),
                clf.predict([[0.5, 0.5]])
            )
            
        except Exception as e:
            self.fail(f"Test failed with error: {str(e)}")
        finally:
            # Cleanup test models directory
            import shutil
            import os
            if os.path.exists("test_models"):
                shutil.rmtree("test_models")

    def test_feature_generation(self):
        # Test feature generator
        sequences = ["ATGCATGC", "GCTAGCTA"]
        generator = feature_generator(sequences)
        features, label = next(generator)
        
        self.assertIsInstance(features, list)
        self.assertIsInstance(label, int)
        self.assertEqual(label, 1)
        
        # Test with invalid sequence
        with self.assertRaises(ValueError):
            next(feature_generator(["ATGN"]))

    def test_incremental_learning(self):
        """Test incremental learning capabilities"""
        try:
            # Create test data chunks
            X1 = np.array([[0.5, 0.5], [0.3, 0.7]])
            y1 = np.array([1, 1])
            X2 = np.array([[0.8, 0.2], [0.1, 0.9]])
            y2 = np.array([0, 0])
            
            # Initialize incremental classifier
            clf = SGDClassifier(random_state=42)
            
            # First chunk
            clf.partial_fit(X1, y1, classes=np.array([0, 1]))
            pred1 = clf.predict(X1)
            self.assertTrue(all(pred1 == y1))
            
            # Second chunk
            clf.partial_fit(X2, y2)
            pred2 = clf.predict(X2)
            self.assertTrue(all(pred2 == y2))
            
            # Test probability predictions
            probs = clf.predict_proba(X1)
            self.assertEqual(probs.shape, (2, 2))  # Two classes
            self.assertTrue(np.all(probs >= 0) and np.all(probs <= 1))
            
        except Exception as e:
            self.fail(f"Test failed with error: {str(e)}")

    def test_uniprot_loading(self):
        try:
            # Test loading UniProt SwissProt database
            sequences, labels = load_uniprot_data(UNIPROT_SPROT_PATH)
            self.assertIsNotNone(sequences)
            self.assertIsNotNone(labels)
            self.assertEqual(len(sequences), len(labels))
            
            # Test sequence validity
            for seq in sequences[:100]:
                self.assertTrue(all(c in 'ATGC' for c in seq))
                self.assertGreaterEqual(len(seq), min(self.VALID_KMER_SIZES))
        
        except Exception as e:
            self.fail(f"UniProt loading failed: {str(e)}")

    def test_supervised_learning_uniprot(self):
        try:
            # Train on UniProt data
            X_train, X_test, y_train, y_test = prepare_uniprot_data(UNIPROT_SPROT_PATH)
            
            # Test model training
            clf = train_kmer_classifier(X_train, y_train)
            self.assertIsInstance(clf, SGDClassifier)
            
            # Evaluate performance
            score = clf.score(X_test, y_test)
            self.assertGreater(score, 0.7)  # Expect at least 70% accuracy
            
            # Test with specific SwissProt sequences
            test_seqs = ["ATGGCCGACTACAAGGACGACGACGACAAG",  # FLAG tag
                        "ATGTACCCATACGATGTTCCAGATTACGCT"]   # HA tag
            for seq in test_seqs:
                pred = predict_kmer(seq, clf)
                self.assertIn(pred, ["Biological", "Artifact"])
                
        except Exception as e:
            self.fail(f"Supervised learning failed: {str(e)}")

    def test_incremental_learning_uniprot(self):
        try:
            # Test incremental learning on UniProt chunks
            data_chunks = get_uniprot_chunks(UNIPROT_SPROT_PATH, chunk_size=1000)
            clf = SGDClassifier(random_state=42)
            
            accuracies = []
            for i, (X_chunk, y_chunk) in enumerate(data_chunks):
                if i == 0:
                    clf.partial_fit(X_chunk, y_chunk, classes=np.array([0, 1]))
                else:
                    clf.partial_fit(X_chunk, y_chunk)
                
                if i % 5 == 0 and i > 0:
                    score = clf.score(X_chunk, y_chunk)
                    accuracies.append(score)
                    
                if i >= 20:  # Test with first 20 chunks
                    break
            
            # Verify learning improvement
            self.assertGreater(accuracies[-1], accuracies[0])
            
        except Exception as e:
            self.fail(f"Incremental learning failed: {str(e)}")

if __name__ == '__main__':
    unittest.main()
