import { createFileRoute } from "@tanstack/react-router";
import { Preprocessor } from "@/components/pre/preprocessor";

export const Route = createFileRoute("/pre")({ component: Preprocessor });
