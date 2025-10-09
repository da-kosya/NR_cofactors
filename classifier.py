# classifier
import json
import os
import re  
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass

@dataclass
class ChainInfo:
    chain_id: str
    description: str
    entity_type: str
    entity_sequence: str

@dataclass
class ClassificationResult:
    pdb_id: str
    is_nr_cofactor_complex: bool
    receptor_chain: Optional[str]
    cofactor_chain: Optional[str]
    receptor_type: Optional[str]
    cofactor_type: Optional[str]
    confidence_score: float
    reasons: List[str]

class NRCofactorClassifier:
    
    def __init__(self):
        # nuclear receptor keywords 
        self.nr_keywords = [
            'nuclear receptor', 'steroid receptor', 'thyroid hormone receptor',
            
            # specific receptors
            'estrogen receptor', 'androgen receptor', 'glucocorticoid receptor',
            'mineralocorticoid receptor', 'progesterone receptor', 
            'peroxisome proliferator-activated receptor', 'ppar',
            'retinoic acid receptor', 'retinoid x receptor', 'rar', 'rxr',
            'liver x receptor', 'lxr', 'farnesoid x receptor', 'fxr',
            'vitamin d receptor', 'vdr', 'thyroid hormone receptor', 'tr',
            'hepatocyte nuclear factor', 'hnf', 'coup-tf', 'nr2f',
            'testicular receptor', 'tr2', 'tr4', 'tailless', 'tlx',
            'photoreceptor cell-specific nuclear receptor', 'pnr',
            'dosage-sensitive sex reversal', 'dax1', 'short heterodimer partner', 'shp'
        ]
        
        # coactivator keywords
        self.coactivator_keywords = [
            # SRC family
            'steroid receptor coactivator', 'src-1', 'src-2', 'src-3', 'src1', 'src2', 'src3',
            'nuclear receptor coactivator', 'ncoa1', 'ncoa2', 'ncoa3', 'grip1', 'tif2', 'actr',
            
            # CBP/p300 family
            'creb-binding protein', 'cbp', 'p300', 'ep300',
            
            # PGC family
            'peroxisome proliferator-activated receptor gamma coactivator', 'pgc-1', 'pgc1',
            'pgc-1alpha', 'pgc-1beta', 'pgc1a', 'pgc1b',
            
            # other 
            'mediator', 'trap', 'drip', 'arc', 'activator',
            'transcriptional intermediary factor', 'receptor-associated protein'
        ]
        
        # corepressor keywords
        self.corepressor_keywords = [
            'nuclear receptor corepressor', 'ncor', 'smrt', 'cornr',
            'silencing mediator', 'corepressor', 'repressor',
            'alien', 'trip15', 'csx', 'hairless', 
        ]

        # self.corepressor_sequence = 'LXXLL' 

    def load_pdb_data(self, pdb_id: str, data_dir: str = "data") -> Optional[Dict]:
        # json files reading
        file_path = os.path.join(data_dir, f"{pdb_id}.json")
        if not os.path.exists(file_path):
            print(f"Warning: Data file not found for {pdb_id}")
            return None
        
        with open(file_path, 'r') as f:
            return json.load(f)
    
    def extract_chain_info(self, data: Dict) -> List[ChainInfo]:
        chains = []
        
        if not data.get('entry') or not data['entry'].get('polymer_entities'):
            return chains
        
        for entity in data['entry']['polymer_entities']:
            # chain IDs
            chain_ids = entity.get('rcsb_polymer_entity_container_identifiers', {}).get('auth_asym_ids', [])
            if not chain_ids:
                continue
            
            # description
            description = entity.get('rcsb_polymer_entity', {}).get('pdbx_description', '').lower()
            
            # entity type
            entity_type = entity.get('entity_poly', {}).get('type', '').lower()
            
            # sequence
            entity_sequence = entity.get('entity_poly', {}).get('pdbx_seq_one_letter_code_can', '')
            
            for chain_id in chain_ids:
                chains.append(ChainInfo(
                    chain_id=chain_id,
                    description=description,
                    entity_type=entity_type,
                    entity_sequence=entity_sequence,
                ))
        
        return chains
    
    def score_nuclear_receptor(self, description: str) -> Tuple[float, List[str]]:
        # how likely the chain is nuclear receptor
        score = 0.0
        reasons = []
        
        desc_lower = description.lower()
        
        for keyword in self.nr_keywords:
            if keyword in desc_lower:
                if keyword in ['nuclear receptor', 'steroid receptor']:
                    score += 1.0
                    reasons.append(f"Generic NR term: {keyword}")
                else:
                    score += 0.8
                    reasons.append(f"Specific NR: {keyword}")
        
        return min(score, 1.0), reasons
    
    def score_cofactor(self, description: str, sequence: str) -> Tuple[float, str, List[str]]:
        # scores how likely a chain is a cofactor based on its description AND sequence.
        coactivator_score = 0.0
        corepressor_score = 0.0
        reasons = []
        
        desc_lower = description.lower()
        
        # coactivator keywords
        for keyword in self.coactivator_keywords:
            if keyword in desc_lower:
                coactivator_score += 0.8
                reasons.append(f"Coactivator keyword: {keyword}")
        
        # corepressor keywords
        for keyword in self.corepressor_keywords:
            if keyword in desc_lower:
                corepressor_score += 0.8
                reasons.append(f"Corepressor keyword: {keyword}")

        
        # search for LXXLL motif in the amino acid sequene
        if sequence and re.search('L..LL', sequence):
            # add a significant score for this strong evidence
            # cannot be small since it is very indicating factor for considering NaR cofactor binding
            coactivator_score += 0.5  
            reasons.append("Corepressor sequence motif 'LXXLL' found")
        
        # determine type and score
        if coactivator_score > corepressor_score:
            return min(coactivator_score, 1.0), "coactivator", reasons
        elif corepressor_score > 0:
            return min(corepressor_score, 1.0), "corepressor", reasons
        else:
            return 0.0, "unknown", reasons
    
    def classify_complex(self, pdb_id: str, data_dir: str = "data") -> ClassificationResult:
            data = self.load_pdb_data(pdb_id, data_dir)
            if not data:
                return ClassificationResult(pdb_id=pdb_id, is_nr_cofactor_complex=False, confidence_score=0.0, reasons=["Data not available"], receptor_chain=None, cofactor_chain=None, receptor_type=None, cofactor_type=None)

            chains = self.extract_chain_info(data)
            if len(chains) < 2:
                return ClassificationResult(pdb_id=pdb_id, is_nr_cofactor_complex=False, confidence_score=0.0, reasons=["Insufficient number of chains"], receptor_chain=None, cofactor_chain=None, receptor_type=None, cofactor_type=None)

            best_pair = {
                "receptor": None,
                "cofactor": None,
                "nr_score": 0.0,
                "cofactor_score": 0.0,
                "confidence": 0.0,
                "nr_reasons": [],
                "cofactor_reasons": []
            }
            
            for i in range(len(chains)):
                for j in range(len(chains)):
                    if i == j:
                        continue

                    receptor_candidate = chains[i]
                    cofactor_candidate = chains[j]
                    
                    nr_score, nr_reasons = self.score_nuclear_receptor(receptor_candidate.description)
                    
                    # pass the sequence to the scoring function
                    cofactor_score, cofactor_type, cofactor_reasons = self.score_cofactor(
                        cofactor_candidate.description, 
                        cofactor_candidate.entity_sequence
                    )
                    
                    penalty_score, _, _ = self.score_cofactor(receptor_candidate.description, receptor_candidate.entity_sequence)
                    nr_score -= penalty_score * 0.5 

                    confidence = (nr_score + cofactor_score) / 2.0

                    if confidence > best_pair["confidence"]:
                        best_pair["receptor"] = receptor_candidate
                        best_pair["cofactor"] = cofactor_candidate
                        best_pair["nr_score"] = nr_score
                        best_pair["cofactor_score"] = cofactor_score
                        best_pair["confidence"] = confidence
                        best_pair["nr_reasons"] = nr_reasons
                        best_pair["cofactor_reasons"] = cofactor_reasons
                        best_pair["cofactor_type"] = cofactor_type

            is_complex = best_pair["nr_score"] >= 0.5 and best_pair["cofactor_score"] >= 0.3
            
            all_reasons = []
            if best_pair["receptor"]:
                all_reasons.extend([f"Receptor ({best_pair['receptor'].chain_id}): " + r for r in best_pair["nr_reasons"]])
            if best_pair["cofactor"]:
                all_reasons.extend([f"Cofactor ({best_pair['cofactor'].chain_id}): " + r for r in best_pair["cofactor_reasons"]])

            if not is_complex and best_pair['confidence'] > 0:
                all_reasons.append(f"Failed threshold check: NR_score({best_pair['nr_score']:.2f}) >= 0.5, Cofactor_score({best_pair['cofactor_score']:.2f}) >= 0.3")

            return ClassificationResult(
                pdb_id=pdb_id,
                is_nr_cofactor_complex=is_complex,
                receptor_chain=best_pair["receptor"].chain_id if best_pair["receptor"] else None,
                cofactor_chain=best_pair["cofactor"].chain_id if best_pair["cofactor"] else None,
                receptor_type=best_pair["receptor"].description if best_pair["receptor"] else None,
                cofactor_type=best_pair.get("cofactor_type") if best_pair["cofactor"] else None,
                confidence_score=best_pair["confidence"],
                reasons=all_reasons
            )
    
    def batch_classify(self, pdb_ids: List[str], data_dir: str = "data") -> List[ClassificationResult]:
        results = []
        for pdb_id in pdb_ids:
            result = self.classify_complex(pdb_id, data_dir)
            results.append(result)
        return results

# Example usage and testing
# def main():
#     classifier = NRCofactorClassifier()
    
#     # Test with your existing data
#     test_ids = ["2QQ1"]  # Add more IDs as you collect data
    
#     results = classifier.batch_classify(test_ids)
    
#     print("\n=== Nuclear Receptor-Cofactor Complex Classification Results ===")
#     print(f"{'PDB ID':<8} {'Complex?':<10} {'Confidence':<12} {'Receptor':<15} {'Cofactor':<15}")
#     print("-" * 80)
    
#     for result in results:
#         complex_status = "YES" if result.is_nr_cofactor_complex else "NO"
#         confidence = f"{result.confidence_score:.2f}"
#         receptor = result.receptor_chain or "N/A"
#         cofactor = result.cofactor_chain or "N/A"
        
#         print(f"{result.pdb_id:<8} {complex_status:<10} {confidence:<12} {receptor:<15} {cofactor:<15}")
        
#         if result.reasons:
#             print("  Reasons:")
#             for reason in result.reasons:
#                 print(f"    - {reason}")
#         print()

# if __name__ == "__main__":
#     main()