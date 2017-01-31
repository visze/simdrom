package de.charite.compbio.simdrom;

import java.io.IOException;
import java.util.Random;

import org.apache.commons.cli.ParseException;

import de.charite.compbio.simdrom.cli.SIMdromSetting;
import de.charite.compbio.simdrom.io.writer.VCFTSVWriter;
import de.charite.compbio.simdrom.sampler.DeNovoSampler;
import de.charite.compbio.simdrom.sampler.SpikeIn;
import de.charite.compbio.simdrom.sampler.vcf.VCFSampleSelector;
import de.charite.compbio.simdrom.sampler.vcf.VCFSampler;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;

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
		VCFSampler.Builder backgroundSamplerBuilder = new VCFSampler.Builder();

		// VCFSampler backgroundSampler = new VCFSampler(SIMdromSetting.BACKGROUND_VCF);
		backgroundSamplerBuilder = backgroundSamplerBuilder.file(SIMdromSetting.BACKGROUND_VCF)
				.probability(SIMdromSetting.BACKGROUND_PROBABILITY);

		// select (random) sample if set
		if (SIMdromSetting.ONLY_ONE_BACKGROUND_SAMPLE) {
			VCFSampleSelector selecter;
			if (SIMdromSetting.ONLY_ONE_BACKGROUND_SAMPLE_NAME == null)
				selecter = new VCFSampleSelector(SIMdromSetting.BACKGROUND_VCF, new Random().nextLong());
			else
				selecter = new VCFSampleSelector(SIMdromSetting.BACKGROUND_VCF,
						SIMdromSetting.ONLY_ONE_BACKGROUND_SAMPLE_NAME);
			backgroundSamplerBuilder.sample(selecter.getSample());
		}
		
		// change the default sample name
		backgroundSamplerBuilder.sampleName(SIMdromSetting.SAMPLE_NAME);
		
		
		
		backgroundSamplerBuilder = backgroundSamplerBuilder
				.afIdentifier(SIMdromSetting.BACKGROUND_ALLELE_FREQUENCY_IDENTIFIER)
				.acIdentifier(SIMdromSetting.BACKGROUND_ALT_ALLELE_COUNT)
				.anIdentifier(SIMdromSetting.BACKGROUND_ALLELE_COUNT);

		backgroundSamplerBuilder = backgroundSamplerBuilder.variantsAmount(SIMdromSetting.BACKGROUND_VARIANT_NUMBER);

		if (SIMdromSetting.INTERVALS.isPresent())
			backgroundSamplerBuilder = backgroundSamplerBuilder.intervals(SIMdromSetting.INTERVALS.get());

		if (SIMdromSetting.USE_DE_NOVO)
			backgroundSamplerBuilder = backgroundSamplerBuilder
					.deNovoGenerator(new DeNovoSampler(SIMdromSetting.DE_NOVO_RATE, SIMdromSetting.REFERENCE));

		// 3) Build writer
		VariantContextWriter writer;
		if (SIMdromSetting.OUTPUT == null)
			writer = new VariantContextWriterBuilder().setOutputVCFStream(System.out)
					.unsetOption(Options.INDEX_ON_THE_FLY).build();
		else
			writer = new VariantContextWriterBuilder().setOutputFile(SIMdromSetting.OUTPUT).build();
		
		// 4) Set VCF for mutation (if set) and settings
		// VCFSampler mutationSampler = null;
		VCFSampler.Builder mutationSamplerBuilder = new VCFSampler.Builder();
		if (SIMdromSetting.MUTATIONS_VCF != null) {

			mutationSamplerBuilder = mutationSamplerBuilder.file(SIMdromSetting.MUTATIONS_VCF)
					.probability(SIMdromSetting.MUTATIONS_PROBABILITY).filters(SIMdromSetting.MUTATIONS_FILTERS);

			// select (random) sample if set
			if (SIMdromSetting.ONLY_ONE_MUTATIONS_SAMPLE) {
				VCFSampleSelector selecter;
				if (SIMdromSetting.ONLY_ONE_MUTATIONS_SAMPLE_NAME == null)
					selecter = new VCFSampleSelector(SIMdromSetting.MUTATIONS_VCF, new Random().nextLong());
				else
					selecter = new VCFSampleSelector(SIMdromSetting.MUTATIONS_VCF,
							SIMdromSetting.ONLY_ONE_MUTATIONS_SAMPLE_NAME);
				mutationSamplerBuilder.sample(selecter.getSample());
			}
			mutationSamplerBuilder = mutationSamplerBuilder
					.afIdentifier(SIMdromSetting.MUTATIONS_ALLELE_FREQUENCY_IDENTIFIER)
					.acIdentifier(SIMdromSetting.MUTATIONS_ALT_ALLELE_COUNT)
					.anIdentifier(SIMdromSetting.MUTATIONS_ALLELE_COUNT);

			mutationSamplerBuilder = mutationSamplerBuilder.variantsAmount(SIMdromSetting.MUTATIONS_VARIANT_NUMBER);

			if (SIMdromSetting.INTERVALS.isPresent())
				mutationSamplerBuilder = mutationSamplerBuilder.intervals(SIMdromSetting.INTERVALS.get());

			// 5) Generate spikein class
			boolean log = SIMdromSetting.SPLIKE_IN_LOGFILE != null;
			SpikeIn spikein = new SpikeIn(backgroundSamplerBuilder.build(), mutationSamplerBuilder.build(), log);
			
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
			
			spikein.close();
		} else {
			VCFSampler sampler = backgroundSamplerBuilder.build();
			writer.writeHeader(sampler.getVCFHeader());
			
			
			while (sampler.hasNext()) {
				VariantContext vc = sampler.next();
				if (vc == null)
					break;
				writer.add(vc);
			}

		}

	

		// 9) close properly and exit properly
		writer.close();
		
		System.exit(0);
	}
}
