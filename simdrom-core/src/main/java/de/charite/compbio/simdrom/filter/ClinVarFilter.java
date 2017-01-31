package de.charite.compbio.simdrom.filter;

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Sets;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.CommonInfo;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

/**
 * An implementation to use the VCF-Info with kex/values provided by the ClinVar file to filter it. It will returns all
 * Variants and alleles that comes from one of the used databases {@link #clndsbn} (OMIM,MedLine,...) that have one of
 * the clinical origins in {@link #clnorigin} (somatic, germline,...) and one of the the clinical significance clinical
 * significances in {@link #clnsig}.
 * 
 * @author <a href="mailto:max.schubach@charite.de">Max Schubach</a>
 *
 */
public class ClinVarFilter extends AFilter {

	/**
	 * All clinical significances that should be used (info field CLNSIG). E.g. 4 and 5 is likely pathogenic and
	 * pathogenic
	 */
	private Set<Integer> clnsig = Sets.newHashSet(4, 5);
	/**
	 * The clinical origin (info field CLNORIGIN). E.g. 1 is germline.
	 */
	private Set<Integer> clnorigin = Sets.newHashSet(1);
	/**
	 * The names in the databases that should be used (info field CLNDSDB). E.g. OMIM
	 */
	private Set<String> clndsbn = Sets.newHashSet("OMIM");

	/**
	 * Default constructor
	 * 
	 * @param clinsig
	 *            Set of clinical significances
	 * @param clinorigin
	 *            Set of clinical origins
	 * @param clndsbn
	 *            Set of clinical databases.
	 */
	public ClinVarFilter(Set<Integer> clinsig, Set<Integer> clinorigin, Set<String> clndsbn) {
		super(FilterType.CLINVAR_FILTER);
		this.clnsig = clinsig;
		this.clnorigin = clinorigin;
		this.clndsbn = clndsbn;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.charite.compbio.simdrom.filter.IFilter#filter(htsjdk.variant. variantcontext.VariantContext)
	 */
	@Override
	public Optional<VariantContext> filter(Optional<VariantContext> optional_vc) {
		if (optional_vc.isPresent()) {
			VariantContext vc = optional_vc.get();
			Map<Integer, HashMultimap<String, String>> sigs_dbs_ids = new HashMap<>();
			List<Allele> newAlleles = new ArrayList<>();
			newAlleles.add(vc.getReference());

			CommonInfo infoField = vc.getCommonInfo();
			if (infoField.hasAttribute("CLNSIG") && infoField.hasAttribute("CLNALLE")
					&& infoField.hasAttribute("CLNDSDB") && infoField.hasAttribute("CLNDSDBID")) {
				Object sigs = infoField.getAttribute("CLNSIG");
				Object alleles = infoField.getAttribute("CLNALLE");
				Object dbs = infoField.getAttribute("CLNDSDB");
				Object dbids = infoField.getAttribute("CLNDSDBID");
				Object origin = infoField.getAttribute("CLNORIGIN");

				if (dbs.getClass().isArray()) {
					filterArray(vc, sigs_dbs_ids, newAlleles, sigs, alleles, dbs, dbids, origin);
				} else if (dbs instanceof String) {

					filterPerClinAllele(vc, sigs_dbs_ids, newAlleles, (String) sigs, (String) dbs, (String) dbids,
							(String) origin, Integer.parseInt((String) alleles));
				} else if (dbs instanceof List)
					filterArray(vc, sigs_dbs_ids, newAlleles, ((List<?>) sigs).toArray(), ((List<?>) alleles).toArray(),
							((List<?>) dbs).toArray(), ((List<?>) dbids).toArray(), ((List<?>) origin).toArray());
				else
					throw new RuntimeException("OHOOOO");
				// no allele matches, return null
				if (newAlleles.size() <= 1)
					return Optional.empty();
				else {
					VariantContextBuilder builder = new VariantContextBuilder().chr(vc.getContig()).start(vc.getStart())
							.noGenotypes().stop(vc.getEnd()).alleles(newAlleles);
					HashMultimap<String, String> dbAttributes = HashMultimap.create();
					for (Integer sig : sigs_dbs_ids.keySet()) {
						for (String db : sigs_dbs_ids.get(sig).keySet()) {
							for (String id : sigs_dbs_ids.get(sig).get(db)) {
								dbAttributes.put(db, id);
							}

						}
					}
					builder = builder.attributes(dbAttributes.asMap());
					builder = builder.attribute("CLNSIG", sigs_dbs_ids.keySet());
					return Optional.of(builder.make());
				}

			}
		}
		return optional_vc;
	}

	private void filterArray(VariantContext vc, Map<Integer, HashMultimap<String, String>> sigs_dbs_ids,
			List<Allele> newAlleles, Object sigs, Object alleles, Object dbs, Object dbids, Object origins) {
		final int dbLength = Array.getLength(dbs);
		for (int i = 0; i < dbLength; i++) { // alleles

			int allelePerCLNAllele = Integer.parseInt((String) Array.get(alleles, i));
			filterPerClinAllele(vc, sigs_dbs_ids, newAlleles, ((String) Array.get(sigs, i)),
					((String) Array.get(dbs, i)), ((String) Array.get(dbids, i)), ((String) Array.get(origins, i)),
					allelePerCLNAllele);

		}
	}

	private void filterPerClinAllele(VariantContext vc, Map<Integer, HashMultimap<String, String>> sigs_dbs_ids,
			List<Allele> newAlleles, String sigs, String dbs, String ids, String origins, int allelePerCLNAllele) {

		String[] sigsPerCLNAllele = sigs.split("\\|");
		String[] dbsPerCLNAllele = dbs.split("\\|");
		String[] idsPerCLNAllele = ids.split("\\|");

		int origin = Integer.parseInt(origins);

		if (allelePerCLNAllele <= 0 || !clnorigin.contains(origin))
			return;

		Allele allele = vc.getAlternateAllele(allelePerCLNAllele - 1);
		boolean use = false;
		for (int j = 0; j < dbsPerCLNAllele.length; j++) { // significance
			int sig = Integer.parseInt(sigsPerCLNAllele[j]);

			if (clnsig.contains(sig)) {
				String[] dbsOfCLNAllele = dbsPerCLNAllele[j].split(":");
				String[] idsOfCLNAllele = idsPerCLNAllele[j].split(":");

				for (int k = 0; k < dbsOfCLNAllele.length; k++) { // databases
					if (clndsbn.contains(dbsOfCLNAllele[k])) {
						if (!sigs_dbs_ids.containsKey(sig)) {
							sigs_dbs_ids.put(sig, HashMultimap.create());
						}
						sigs_dbs_ids.get(sig).put(dbsOfCLNAllele[k], idsOfCLNAllele[k]);
						use = true;
					}
				}

			}

		}
		if (use) {
			if (!newAlleles.contains(allele))
				newAlleles.add(allele);
		}
	}

}
