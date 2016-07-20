package de.charite.compbio.simdrom;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;

import java.io.IOException;

import org.apache.commons.cli.ParseException;

import de.charite.compbio.simdrom.cli.SIMdromSetting;
import de.charite.compbio.simdrom.io.writer.VCFTSVWriter;
import de.charite.compbio.simdrom.sampler.DeNovoSampler;
import de.charite.compbio.simdrom.sampler.SpikeIn;
import de.charite.compbio.simdrom.sampler.vcf.VCFRandomSampleSelecter;
import de.charite.compbio.simdrom.sampler.vcf.VCFSampler;

/**
 * Main class for the command line interface.
 * 
 * @author <a href="mailto:max.schubach@charite.de">Max Schubach</a>
 *
 */
public class Main {

	public static void main(String[] args) throws ParseException, IOException {

		// 1) Parse options
		SIMdromSetting.parse(args);

		// 2) Set VCF for background population and settings
		VCFSampler backgroundSampler = new VCFSampler(SIMdromSetting.BACKGROUND_VCF);

		backgroundSampler.setProbability(SIMdromSetting.BACKGROUND_PROBABILITY);
		if (SIMdromSetting.ONLY_ONE_SAMPLE) {
			VCFRandomSampleSelecter selecter;
			if (SIMdromSetting.ONLY_ONE_SAMPLE_NAME == null)
				selecter = new VCFRandomSampleSelecter(SIMdromSetting.BACKGROUND_VCF);
			else 
				selecter = new VCFRandomSampleSelecter(SIMdromSetting.BACKGROUND_VCF, SIMdromSetting.ONLY_ONE_SAMPLE_NAME);
			backgroundSampler.setSample(selecter.getSample());
		}
		if (SIMdromSetting.BACKGROUND_ALLELE_FREQUENCY_IDENTIFIER != null) {
			backgroundSampler.setAFIdentifier(SIMdromSetting.BACKGROUND_ALLELE_FREQUENCY_IDENTIFIER);
		}
		if (SIMdromSetting.BACKGROUND_ALT_ALLELE_COUNT != null && SIMdromSetting.BACKGROUND_ALLELE_COUNT != null) {
			backgroundSampler.setACIdentifier(SIMdromSetting.BACKGROUND_ALT_ALLELE_COUNT);
			backgroundSampler.setANIdentifier(SIMdromSetting.BACKGROUND_ALLELE_COUNT);
		}
		if (SIMdromSetting.BACKGROUND_VARIANT_NUMBER > 0) {
			backgroundSampler.setVariantsAmount(SIMdromSetting.BACKGROUND_VARIANT_NUMBER);
		}
		if (SIMdromSetting.INTERVALS != null)
			backgroundSampler.setIntervals(SIMdromSetting.INTERVALS);
		if (SIMdromSetting.USE_DE_NOVO)
			backgroundSampler.setDeNovoGenerator(new DeNovoSampler(SIMdromSetting.DE_NOVO_RATE,SIMdromSetting.REFERENCE));

		// 3) Set VCF for mutation (if set) and settings
		VCFSampler mutationSampler = null;
		if (SIMdromSetting.MUTATIONS_VCF != null) {
			mutationSampler = new VCFSampler(SIMdromSetting.MUTATIONS_VCF);
			mutationSampler.setFilters(SIMdromSetting.MUTATIONS_FILTERS);
			mutationSampler.setProbability(SIMdromSetting.MUTATIONS_PROBABILITY);
			if (SIMdromSetting.MUTATIONS_ALLELE_FREQUENCY_IDENTIFIER != null) {
				mutationSampler.setAFIdentifier(SIMdromSetting.MUTATIONS_ALLELE_FREQUENCY_IDENTIFIER);
			}
			if (SIMdromSetting.MUTATIONS_VARIANT_NUMBER > 0) {
				mutationSampler.setVariantsAmount(SIMdromSetting.MUTATIONS_VARIANT_NUMBER);
			}
			if (SIMdromSetting.MUTATIONS_ALT_ALLELE_COUNT != null && SIMdromSetting.MUTATIONS_ALLELE_COUNT != null) {
				mutationSampler.setACIdentifier(SIMdromSetting.MUTATIONS_ALT_ALLELE_COUNT);
				mutationSampler.setANIdentifier(SIMdromSetting.MUTATIONS_ALLELE_COUNT);
			}
			if (SIMdromSetting.INTERVALS != null)
				mutationSampler.setIntervals(SIMdromSetting.INTERVALS);
		}

		// 4) Build writer
		VariantContextWriter writer;
		if (SIMdromSetting.OUTPUT == null) 
			writer = new VariantContextWriterBuilder().setOutputVCFStream(System.out)
				.unsetOption(Options.INDEX_ON_THE_FLY).build();
		else
			writer = new VariantContextWriterBuilder().setOutputFile(SIMdromSetting.OUTPUT).build();

		// 5) Generate spikein class
		boolean log = SIMdromSetting.SPLIKE_IN_LOGFILE != null;
		SpikeIn spikein = new SpikeIn(backgroundSampler, mutationSampler, log);

		// 6) write out VCF header
		writer.writeHeader(spikein.getVCFHeader());

		// 7) spike in and write out
		while (spikein.hasNext()) {
			VariantContext vc = spikein.next();
			if (vc == null)
				break;
			writer.add(vc);
		}

		// 8) write log if set
		if (log) {
			VCFTSVWriter logWriter = new VCFTSVWriter(SIMdromSetting.SPLIKE_IN_LOGFILE);
			boolean header = false;
			for (VariantContext vc : spikein.getVcLogs()) {
				if (!header) {
					logWriter.writeHeader(vc);
					header = true;
				}
				logWriter.add(vc);
			}
			logWriter.close();
		}

		// 9) close properly and exit properly
		writer.close();
		spikein.close();
		System.exit(0);
	}
}
